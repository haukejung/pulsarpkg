import os
import sqlite3
import glob
import numpy as np
import warnings

class Files:
    """
    Files class for accessing fits files directly
    """
    def import_fits(self):
        if not self.imported:
            import astropy.io.fits as fits
            self.fits = fits
            self.imported = True

    def __init__(self, file, debug=False, verbose=False):
        self.file = file[0]
        self.debug = debug
        self.imported = False
        self.fits = None
        self.files = self.get_file_list(file)
        self.verbose = verbose
        if not debug:
            warnings.simplefilter('ignore', UserWarning)

    def get_file_list(self, search_list, debug=False, verbose=False):
        """
        Get the file list based on the search_list.
        :param search_list: list of filenames, can include wildcards, e.g. dir/*.fits
        :return: List of files
        """
        # self.import_fits()
        self.import_fits()
        file_list = []
        for string in search_list:
            for stri in glob.glob(string):
                file_list.append(stri)
        # assert file_list != [], 'No files found.'
        return file_list

    def get_header_data(self, filename):
        """
        Get the header information and data of the provided file.
        :param filename:
        """
        self.import_fits()
        hdulist = self.fits.open(filename)
        hdulist.verify('fix')
        header = hdulist[0].header
        astrodata = hdulist[0].data

        header = self.fix_header(header)
        hdulist[0].header = header
        hdulist[0].data = astrodata
        # hdulist.close()
        return hdulist, header, astrodata

    def fix_header(self, header):
        """
        Fix header fields: remove empty headers, set None where no familiar value is specified
        :param header:
        :return:
        """
        self.import_fits()
        if '' in header.keys():
            del(header[''])     # removing all empty header
        for key, value in header.items():
            if isinstance(value, self.fits.card.Undefined):
                header[key] = None
                if self.debug:
                    print('fixed', key)
        # del(header[""])     # removing all empty headers
        return header


class DB(Files):
    """
    FitsDB class for storing fits files in a sqlite database
    """
    def __init__(self, file, debug=False, verbose=False):
        assert isinstance(file, list)
        Files.__init__(self, file, debug, verbose)
        self.fraction = 0
        db = os.access(self.file, os.F_OK)
        self.conn = sqlite3.connect(self.file, detect_types=sqlite3.PARSE_DECLTYPES)
        self.conn.row_factory = sqlite3.Row  # makes the results of querys a dict instead of a tuple
        self.cursor = self.conn.cursor()
        if os.path.isfile(self.file):
            if not db:
                self.create_table()

        self.cursor.execute("PRAGMA locking_mode=EXCLUSIVE;")
        sqlite3.enable_callback_tracebacks(True)

    def __del__(self):
        if self.conn:
            self.conn.close()

    def create_table(self):
        """
        Create the SQLite tables
        """
        command = 'CREATE TABLE IF NOT EXISTS headers (' \
                  'id INTEGER PRIMARY KEY, ' \
                  'filename TEXT UNIQUE ON CONFLICT REPLACE, ' \
                  'ctime DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP, ' \
                  'mtime DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP, ' \
                  'keywords TEXT DEFAULT NULL);'
        self.cursor.execute(command)

        command = 'CREATE TRIGGER IF NOT EXISTS headers_update_trigger AFTER UPDATE ON headers FOR EACH ROW ' \
                  'BEGIN ' \
                  'UPDATE headers SET mtime=CURRENT_TIMESTAMP WHERE id=NEW.id; ' \
                  'END'
        self.cursor.execute(command)

        command = 'CREATE TABLE IF NOT EXISTS astrodata (' \
                  'headers_id INTEGER REFERENCES headers(id) ON DELETE CASCADE, '\
                  'DATA BLOB);'
        self.cursor.execute(command)

        self.conn.commit()

    def get_id(self, filename):
        """
        Return the id of the row that matches that filename. Return [] if there is no match. Raise AssertionError if
        there is more than one match.
        :param filename: Search for the row with the given filename
        :return: id of the matching row, or []
        """
        if self.debug:
            print('filename', filename)
        command = 'SELECT id FROM headers WHERE filename = "{0}"'.format(os.path.basename(filename))
        self.cursor.execute(command)
        header_id = self.cursor.fetchall()
        assert len(header_id) in [0, 1], 'Multiple rows match filename.'
        return header_id

    def sql(self, command='SELECT * FROM headers;'):
        """
        Issue an arbitrary SQL-Statement
        :param command: SQL-statement
        :return: matching rows
        """
        self.cursor.execute(command)
        res = self.cursor.fetchall()
        return res

    def get_columns(self):
        """
        Get a list of columns from the headers table
        :return: List of columns
        """
        self.cursor.execute('PRAGMA table_info(headers)')
        res = self.cursor.fetchall()
        return [re[1] for re in res]  # re[0-4] are column number, name, type, notnull, default value, PRIMARY_KEY

    def check_columns(self, headers):
        """
        Checks the database if the header exists as a column in the database
        If not, it adds a column
        :param headers: list of headers
        :rtype: None
        """
        self.import_fits()
        columns = self.get_columns()
        for header in headers:
            if header not in columns and header is not '':
                self.cursor.execute('ALTER TABLE headers ADD "%s" NUMERIC' % header)
                columns.append(header)
        self.conn.commit()

    def report_percentage(self):
        return round(self.fraction*100)

    def ingest_data(self, search_list):
        """
        Ingest the header information and data with an UPDATE or INSERT.
        :param search_list:
        """
        self.import_fits()
        files = self.get_file_list(search_list)
        if self.debug:
            print(files)
        self.fraction = 0
        for file in files:
            header_id = self.get_id(file)     # will be [] if file is not found in the DB
            hdulist, header, astrodata = self.get_header_data(file)
            # del(header[""])     # removing all empty headers
            # print('header', header)
            self.check_columns(header.keys())
            if not header_id:
                # INSERT
                keys = []
                keys.extend(['filename', 'keywords'])
                keys.extend(header.keys())

                command = 'INSERT INTO headers (`'
                command += '`, `'.join(keys) + '`)'
                # command += 'filename, ctime, mtime, keywords'
                # for key in header.keys():
                #     command += '"' + key + '",'

                values = []
                values.extend([os.path.basename(file), os.path.relpath(file)])
                # values.extend([v for v in header.values()])
                values.extend(header.values())

                command += ' VALUES (' + \
                           '?, '*(len(header)-1+2) + '?)'        # +2 for the filename and keywords fields

                if self.debug:
                    log_sql_stmt(self.cursor, command, *values)
                self.cursor.execute(command, tuple(values))
                self.conn.commit()
                # command = 'INSERT INTO astrodata (headers_id, dim1, dim2, DATA) VALUES (?, ?, ?, ?)'
                # self.cursor.execute(command, (self.cursor.lastrowid, len(astrodata), len(astrodata[0]), astrodata))
                command = 'INSERT INTO astrodata (headers_id, DATA) VALUES (?, ?)'
                self.cursor.execute(command, (self.cursor.lastrowid, astrodata))
                self.conn.commit()
            else:
                # UPDATE
                command = 'UPDATE OR ROLLBACK headers SET '
                for key in header:
                    command += '"' + key + '" = "' + str(header[key]) + '",'
                command = command[:-1]
                command += ' WHERE id = "' + str(header_id[0][0]) + '"'     # get first row and first field (id)
                self.cursor.execute(command)
                self.conn.commit()
                # command = 'UPDATE astrodata SET dim1 = ?, dim2 = ?, DATA = ? WHERE headers_id = ?'
                # self.cursor.execute(command, (astrodata, len(astrodata), len(astrodata[0]), str(header_id[0][0])))
                command = 'UPDATE astrodata SET DATA = ? WHERE headers_id = ?'
                self.cursor.execute(command, (astrodata, str(header_id[0][0])))
                self.conn.commit()
            self.fraction += 1/len(files)
            if self.verbose:
                print('\r{0}%'.format(self.report_percentage()), end='')
        print("")

    def extract(self, attributes: dict):
        """
        Extract data from the database
        :type attributes: dict
        :return: List of astropy HDU lists
        """
        self.import_fits()
        columns = self.get_columns()
        command = 'SELECT * FROM headers '
        command += 'JOIN astrodata ON headers.id = astrodata.headers_id'  # also get the dynamic spectrum

        values = []
        if len(attributes) > 0:  # if we have attributes to search for, add the WHERE clause and the attributes
            command += ' WHERE '
            for attribute, value in attributes.items():
                # check if attributes are present in the database
                if attribute not in columns:
                    raise NameError("Attribute %s not found in database" % attribute)
                if attribute is 'MJD':
                    command += '`MJD` >= ? AND `MJD` <= ? AND '
                    values.extend(value.split(" "))
                elif attribute is 'FREQ':
                    command += '`FREQ` >= ? AND `FREQ` <= ? AND '
                    values.extend(value.split(" "))
                else:
                    split = value.split(" ")
                    if len(split) == 2:
                        command += '`{0}` >= ? AND `{1}` <= ? AND '.format(attribute, attribute)
                        values.extend(value.split(" "))
                    elif len(split) == 1:
                        command += '{0} LIKE ? AND '.format(attribute)
                        values.extend(['%'+value+'%'])
                    else:
                        raise ValueError("Didn't understand values for attribute {0}".format(attribute))
            command = command[:-4]  # remove the last "AND "

        if self.debug:
            log_sql_stmt(self.cursor, command, *values)
        self.cursor.execute(command, values)
        rows = self.cursor.fetchall()

        ret_list = []
        for row in rows:
            hdulist = self.create_hdulist_from_row(row)
            # TODO: could delete the DATA field as we don't need it from now on
            ret_list.append((hdulist, row))
        return ret_list

    def create_hdulist_from_row(self, row, writetofile=False):
        """
        Creates an astropy HDUList object from a database row
        :param row: database row
        :return: HDUList object
        """
        header = self.fits.Header()
        for key in row.keys():
            key = str(key)
            if key in ['id', 'filename', 'ctime', 'mtime', 'keywords', 'headers_id', 'DATA']:
                # we don't want those in the header
                continue
            header.extend([(key, row[key])])
        hdulist = self.fits.HDUList()  # start creating the new HDU list

        dtype = np.float32 if row["BITPIX"] == -32 else np.float64  # BITPIX shows if floats are 32 or 64 bit long
        data = np.fromstring(row["DATA"], dtype=dtype)
        if row["BITPIX"] == -32:
            data = data.byteswap()  # change endianness only if 32-bit floats
        if len(data) == row["NAXIS1"]*row["NAXIS2"]:
            data = data.reshape((row["NAXIS2"], row["NAXIS1"]))  # recreate the 2D matrix from the 1D array
        else:
            print("The NAXIS parameters don't match the data (%i <> %i)."
                  % (len(data), row["NAXIS2"]*row["NAXIS1"]))
            return hdulist  # return the empty list to not abort the program
        data = np.rot90(data)

        imagehdu = self.fits.ImageHDU(data=data, header=header)
        hdulist.append(imagehdu)

        if writetofile:
            hdulist.writeto(row["filename"], output_verify='fix')  # writes the file back to a file
        return hdulist


def log_sql_stmt(cursor, sql, *args):
    if len(args) > 0:
        # generates SELECT quote(?), quote(?), ...
        cursor.execute("SELECT " + ", ".join(["quote(?)"]*len(args)), args)
        quoted_values = cursor.fetchone()
        for quoted_value in quoted_values:
            sql = sql.replace('?', quoted_value, 1)
    print("SQL command: " + sql)

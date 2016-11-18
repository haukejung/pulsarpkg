from __future__ import print_function

import functions

import os
import sqlite3
import glob
import numpy as np
import warnings
import io


def adapt_array(arr):
    """
    Converts an array to bytes
    :param arr: The array to convert
    :return: The sqlite3 binary
    """
    out = io.BytesIO()
    np.save(out, arr)
    out.seek(0)
    return sqlite3.Binary(out.read())


def convert_array(text):
    """
    Convert bytes to numpy array
    :param text: binary text
    :return: numpy array
    """
    out = io.BytesIO(text)
    out.seek(0)
    return np.load(out)


class Files:
    """
    Files class for accessing fits files directly
    """
    def import_fits(self):
        if not self.imported_fits:
            try:
                import astropy.io.fits as fits
                self.fits = fits
                self.imported_fits = True
            except ImportError:
                print('For this type of action, the python package "astropy" is needed.')
                exit(1)

    def import_psrchive(self):
        if not self.imported_psrchive:
            try:
                import psrchive
                self.psrchive = psrchive
                self.imported_psrchive = True
            except ImportError:
                print('For this type of action, psrchive must be installed, and the python interface available.\n'
                      'PSRchive also only works with python2')
                exit(1)

    def __init__(self, file, debug=False, verbose=False):
        self.file = file[0]
        self.debug = debug
        self.imported_fits = False
        self.imported_psrchive = False
        self.fits = None
        self.psrchive = None
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

    def get_data_fits(self, filename):
        """
        Get the header information and data of the provided fits-file.
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

    def get_data_ar(self, filename):
        """
        Get the header information and data of the provided psrchive-file. Converts the data to a dynamic spectrum
        :param filename:
        """
        self.import_fits()
        self.import_psrchive()
        ar = self.psrchive.Archive_load(filename)
        ar.pscrunch()
        dedispersed = ar.get_dedispersed()
        if not dedispersed:  # dedisperse the data if it's not already
            ar.dedisperse()

        ar.remove_baseline()
        ar.get_filename()

        # Get metadata
        backend = ar.get_backend_name()
        freq = ar.get_centre_frequency()
        bw = ar.get_bandwidth()
        nbin = ar.get_nbin()
        nchan = ar.get_nchan()
        nsubint = ar.get_nsubint()
        int_len = ar.integration_length()
        source = ar.get_source()
        origin = ar.get_telescope()
        mjd = ar.get_Integration(0).get_start_time().in_days()
        dm = ar.get_dispersion_measure()

        # Using the ProfileShiftFit class to compute the SNR for every subint/channel
        # (analogous to the dynamic_spectra.C code)
        prof_shift = self.psrchive.ProfileShiftFit()
        prof_shift.choose_maximum_harmonic = True
        tot = ar.total()
        tot_prof = tot.get_Profile(0, 0, 0)
        prof_shift.set_standard(tot_prof)

        dyn = np.empty(shape=[nsubint, nchan])
        for i in range(nsubint):
            for j in range(nchan):
                profile = ar.get_Profile(i, 0, j)
                prof_shift.set_Profile(profile)
                dyn[i, j] = prof_shift.get_snr()*profile.get_weight()

        header = {'FREQ': freq, 'BW': bw, 'NCHAN': nchan, 'NSUB': nsubint, 'T_INT': int_len, 'SOURCE': source,
                  'ORIGIN': origin, 'MJD': mjd, 'NAXIS': 2, 'NAXIS1': len(dyn[0]), 'NAXIS2': len(dyn), 'BITPIX': 32, 'DM': dm}

        fits_header = self.fits.Header()  # prepare FITS header
        for key in header:
            # print(key)
            fits_header.extend([(str(key), header[key])])

        hdulist = self.fits.HDUList()  # start creating the new HDU list
        imagehdu = self.fits.ImageHDU(data=dyn, header=fits_header)
        hdulist.append(imagehdu)

        hdulist.writeto(filename[:filename.rfind('.')]+".fits", output_verify='fix')  # writes the file back to a file
        # and changes existing extension to .fits

        return hdulist, header, dyn

    def get_data(self, file):
        ext = os.path.splitext(file)[1]
        hdulist, header, astrodata = None, None, None
        if ext in ['.fit', '.fits']:
            try:
                hdulist, header, astrodata = self.get_data_fits(file)
            except OSError:
                print('{0} is not a FITS file'.format(file))
                exit(1)
        elif ext in ['.ar']:
            hdulist, header, astrodata = self.get_data_ar(file)
            # try:
            #     hdulist, header, astrodata = self.get_data_ar(file)
            # except Exception:
            #     print('File is not in psrchive format')
            #     exit(1)
        else:
            print('Testing if it\'s a pulsar archive..')
            hdulist, header, astrodata = self.get_data_ar(file)
            # try:
            #     print('Testing if it\'s a pulsar archive..')
            #     hdulist, header, astrodata = self.get_data_ar(file)
            # except Exception:
            #     print('File is neither in FITS nor in psrchive format')
            #     exit(1)
        # del(header[""])     # removing all empty headers
        # print('header', header)
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
        functions.check_object_type(file, list)
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

            hdulist, header, astrodata = self.get_data(file)

            self.check_columns(header.keys())
            if not header_id:  # file not in the database yet
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
                self.cursor.execute(command, (self.cursor.lastrowid, adapt_array(astrodata)))
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
                self.cursor.execute(command, (adapt_array(astrodata), str(header_id[0][0])))
                self.conn.commit()
            self.fraction += 1/len(files)
            if self.verbose:
                print('\r{0}%'.format(self.report_percentage()), end='')

    def extract(self, attributes, writetofile=False):
        """
        Extract data from the database
        :type attributes: dict
        :return: List of astropy HDU lists
        """
        functions.check_object_type(attributes, dict)

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
                    raise NameError("Attribute {0} not found in database".format(attribute))
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
            hdulist = self.create_hdulist_from_row(row, writetofile)
            # TODO: could delete the DATA field as we don't need it from now on
            ret_list.append((hdulist, row))
        return ret_list

    def delete(self, id_list):
        if len(id_list) == 0:
            raise ValueError('No entries to delete')
        command = 'DELETE FROM astrodata WHERE headers_id = ?'
        command2 = 'DELETE FROM headers WHERE id = ?'
        for id_ in range(1, len(id_list)):
            command += ' or headers_id = ?'
            command2 += ' or id = ?'

        if self.debug:
            log_sql_stmt(self.cursor, command, *id_list)
            log_sql_stmt(self.cursor, command2, *id_list)
        self.cursor.execute(command, tuple(id_list))
        self.cursor.execute(command2, tuple(id_list))
        self.conn.commit()

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
                # we don't want those in the FITS-header
                continue
            header.extend([(key, row[key])])
        hdulist = self.fits.HDUList()  # start creating the new HDU list

        data = convert_array(row["DATA"])  # , dtype=dtype)
        # dtype = np.float32 if row["BITPIX"] == -32 else np.float64  # BITPIX shows if floats are 32 or 64 bit long
        # if row["BITPIX"] == -32:
        #     data = data.byteswap()  # change endianness only if 32-bit floats
        # if len(data) == row["NAXIS1"]*row["NAXIS2"]:
        #     data = data.reshape((row["NAXIS2"], row["NAXIS1"]))  # recreate the 2D matrix from the 1D array
        if len(data)*len(data[0]) != row["NAXIS1"]*row["NAXIS2"]:
            print("The NAXIS parameters don't match the data ({0}, {1} <> {2}, {3}).".format(
                    len(data), len(data[0]), row["NAXIS2"]*row["NAXIS1"]))
            return hdulist  # return the empty list to not abort the program

        imagehdu = self.fits.ImageHDU(data=data, header=header)
        hdulist.append(imagehdu)
        if writetofile:
            filename = row["filename"]
            ext = os.path.splitext(filename)[1]
            if ext not in ['.fit', '.fits']:
                filename += '.fits'
            hdulist.writeto(filename, output_verify='fix')  # writes the hdulist back to a FITS-file

        return hdulist


def log_sql_stmt(cursor, sql, *args):
    if len(args) > 0:
        # generates SELECT quote(?), quote(?), ...
        cursor.execute("SELECT " + ", ".join(["quote(?)"]*len(args)), args)
        quoted_values = cursor.fetchone()
        for quoted_value in quoted_values:
            sql = sql.replace('?', quoted_value, 1)
    print("SQL command: " + sql)

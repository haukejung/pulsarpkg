from __future__ import division
from __future__ import with_statement
from __future__ import print_function

from fitsdb import sqlite
from arcfinder import computing
from arcfinder import plotting

import astropy.utils.console

import argparse
from collections import OrderedDict


def parse_args():
    """
    Parse the command line arguments.
    """
    parser = argparse.ArgumentParser(
            description='Inserts FITS header information in a SQLite database.')
    subparsers = parser.add_subparsers()
    ingest = subparsers.add_parser('ingest', help='ingest files')
    ingest.set_defaults(subcmd='ingest')
    ingest.add_argument('files', help='filenames, e.g. "dir/*.fits"', nargs='+')
    sql = subparsers.add_parser('sql', help='SQL-Statement to query, e.g. "SELECT * FROM headers")')
    sql.set_defaults(subcmd='sql')
    query = subparsers.add_parser('query', help='Query the database')
    query.set_defaults(subcmd='query')
    query.add_argument('--filename', help='filename')
    query.add_argument('--pulsar', help='Full pulsar name')
    query.add_argument('--mjd', help='MJD-range in the form of "51000 52000"')
    query.add_argument('--freq', help='Frequency bandwidth (MHz) in the form of "300 400"')
    query.add_argument('-a', '--attr', help='other attribute, e.g. --attr "T_INT"', nargs='+')
    query.add_argument('-r', '--attr-value', help='attribute value or range (see --mjd)', nargs='+')
    query.add_argument('--attr-list', action="store_true", help='get a list of available attributes')
    query.add_argument('-d', '--dyn', action="store_true", help='Plot the dynamic spectrum')
    query.add_argument('-s', '--sec', action="store_true", help='Plot the secondary spectrum')
    query.add_argument('-t', '--store', action="store_true", help='store the images to separate files,\n'
                                                                  'defaults to .png if --format is not specified')
    query.add_argument('-p', '--pdf', action="store_true", help='store all of the images in a single pdf file')
    query.add_argument('-m', '--format', help='format for image saving, i.e. "jpg", "eps"')
    query.add_argument('-wf', '--write-files', action="store_true", help='write the rows from the DB back to the files')
    query.add_argument('--csv', action="store_true", help='Write a csv file')
    query.add_argument('--cmap', help='Choose the colormap from the matplotlib palette, default is "viridis"')
    query.add_argument('-e', '--delete', action="store_true", help='Delete the entries that match the query\n'
                                                                   'DOESN\'T ASK FOR CONFIRMATION')

    positional = parser.add_mutually_exclusive_group(required='True')
    positional.add_argument('-b', '--db', action="store_true", help='switch for using a sqlite database')
    positional.add_argument('-f', action="store_true", help='switch for using local files')
    parser.add_argument('file', help='path to the file(s)', nargs='+')
    parser.add_argument('-v', '--verbose', action="store_true", help='enable verbose mode')
    parser.add_argument('--debug', action="store_true", help='enable debug mode')
    args = parser.parse_args()
    if args.debug:
            print(args)
    return args


def create_attr_dict(args, outp, stdlength=10):
    """
    Creates the attr_dict for the db.extract() function
    adds output attributes if additional were specified on the commandline
    :param args: cmdline arguments
    :param outp: attributes to output and corresponding text widths
    :return: attr_dict
    """
    attr_dict = {}
    if args.filename:
        attr_dict["filename"] = args.filename
    if args.pulsar:
        attr_dict["SOURCE"] = args.pulsar
    if args.mjd:
        attr_dict["MJD"] = args.mjd
    if args.freq:
        attr_dict["FREQ"] = args.freq

    # add output attributes if specified
    if args.attr and args.attr_value:
        for a, v in zip(args.attr, args.attr_value):
            a = a.strip('\'\"')
            v = v.strip('\'\"')
            attr_dict[a] = v
            if a not in outp.keys():
                outp[a] = stdlength  # std length for attributes is 10
    return attr_dict


def main(args):
    """
    The main controller.
    :param args: arguments provided by the commandline
    """

    db = None
    if args.db:
        db = sqlite.DB(args.file, args.debug, args.verbose)
    elif args.f:
        files = sqlite.Files(args.file, args.debug, args.verbose)
    if args.subcmd == 'ingest':
        file_list = db.get_file_list(args.files)
        with astropy.utils.console.ProgressBarOrSpinner(len(file_list), "Ingesting..") as bar:
            for i in range(1, len(file_list)):
                db.ingest_data([file_list[i]])
                bar.update(i)
    elif args.subcmd == 'sql':
        for row in db.sql(args.sql):
            print(tuple(row))
    elif args.subcmd == 'query':
        if args.cmap:
            plotting.cmap = args.cmap

        if args.attr_list:
            attrs = db.get_columns()
            print("List of available attributes:")
            for attr in attrs:
                print(attr)
            exit(0)

        # outp = OrderedDict(id=4, filename=30, ORIGIN=23, MJD=18, FREQ=12, BW=10)
        outp = OrderedDict()  # consists of the output columns and their corresponding widths
        if not args.f:
            outp['id'] = 4
        outp['filename'] = 30
        outp['ORIGIN'] = 23
        outp['MJD'] = 18
        outp['FREQ'] = 12
        outp['BW'] = 10

        attr_dict = create_attr_dict(args, outp)
        if args.debug:
            print('attr_dict', attr_dict)

        if args.db:
            result = db.extract(attr_dict, args.write_files)
        elif args.f:
            result = files.files
            pass
        else:
            result = []
        type = 'rows' if args.db else 'files'
        print("Found {0} matching {1}\n".format(len(result), type))

        if len(result) < 1:
            exit(0)  # nothing to do
        else:  # create list header for console output
            text, hline = '', ''
            for key, value in outp.items():
                text += ' {0:{1}} |'.format(key, value)
                hline += '{0}|'.format('-'*(value+2))  # makes the hline under the header; +2 for the two spaces
            print(text[:-2])  # remove the last "\t|"
            print(hline[:-1])  # remove the last "|"

        csv = ''

        all_keys = []  # gather headers of files for csv in there
        if args.f:
            all_keys.append("filename")
        pdf = None
        if args.pdf:
            pdf = plotting.Pdf(attr_dict)

        delete_ids = []
        for res in result:
            if args.db:
                if args.delete:
                    delete_ids.append(res[1]["id"])
                rotate = False
                hdulist = res[0]
                header = res[1]
                filename = header['filename']
                rotate = True
            elif args.f:
                hdulist, header, data = files.get_header_data(res)
                filename = res
                rotate = True if not args.write_files else False

            # csv
            if args.csv:
                if args.db:
                    for key in header.keys():
                        if key != 'DATA':
                            csv += '{0};'.format(header[key])
                    csv += '\n'
                elif args.f:
                    csv += '{0};'.format(res)
                    for key in header.keys():
                        if key not in all_keys:
                            all_keys.append(key)
                            csv += '{0};'.format(header[key]) if key in all_keys else ';'
                    csv += '\n'

            # text output
            text = ''
            if args.f:
                text += ' {0:{1}} |'.format(filename, outp['filename'])
            for key, value in outp.items():
                if args.f and key == 'filename':
                    continue
                text += ' {0:{1}} |'.format(header[key], value)
            print(text[:-2])  # remove the last " |"

            # Plotting:
            if args.dyn and not args.sec:  # only plot dynamic
                dyn = computing.Dynamic(hdulist, header, filename, rotate)
                plotting.show_dyn(dyn, args.store, args.format, pdf)
            if args.sec:  # plot dynamic and secondary
                sec = computing.Secondary(hdulist, header, filename, rotate)
                if args.dyn:
                    plotting.show_dyn(sec, args.store, args.format, pdf)
                plotting.show_sec(sec, args.store, args.format, pdf)
            if not args.store:  # don't plot to screen when storing images
                plotting.show()
            # end for-loop

        if args.delete:
            db.delete(delete_ids)  # delete all ids in "delete"
            print('\nDeleted all matching rows!')

        if args.csv:
            csv_head = ''
            if args.db:
                for key in result[0][1].keys():
                    if key != 'DATA':
                        csv_head += '{0};'.format(key)
                csv_head += '\n'
            if args.f:
                for key in all_keys:
                    csv_head += '{0};'.format(key)
                csv_head += '\n'
            name = 'pulsarpkg_query'  # filename
            for attr, value in attr_dict.items():
                name += '_{0}_{1}'.format(attr, value.replace(' ', '_'))
            name += '_{0}.csv'.format(plotting.datenow)
            with open(name, 'w') as f:
                f.write(csv_head+csv)


if __name__ == '__main__':
    args = parse_args()
    main(args)



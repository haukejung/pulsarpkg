#! /usr/bin/env python
"""

Command line help available by running:
> python fits2sqlite -h

"""

import argparse
from fitsdb import sqlite
from arcfinder import computing
from arcfinder import plotting


def parse_args():
    """
    Parse the command line arguments.
    """
    parser = argparse.ArgumentParser(
            description='Inserts FITS header information in a SQLite database.')
    subparsers = parser.add_subparsers()
    ingest = subparsers.add_parser('ingest', help='ingest files')
    ingest.set_defaults(subcmd='ingest')
    ingest.add_argument('files', help='filenames, e.g. "dir/*.fits"')
    sql = subparsers.add_parser('sql', help='SQL-Statement to query, e.g. "SELECT * FROM headers")')
    sql.set_defaults(subcmd='sql')
    query = subparsers.add_parser('query', help='Query the database')
    query.set_defaults(subcmd='query')
    query.add_argument('--filename', help='filename')
    query.add_argument('--pulsar', help='Full pulsar name')
    query.add_argument('--mjd', help='MJD-range in the form of "51000 52000"')
    query.add_argument('--freq', help='Frequency bandwidth (MHz) in the form of "300 400"')
    query.add_argument('-a-', '--attr', help='other attribute, e.g. --attr "T_INT', nargs='*')
    query.add_argument('-av', '--attr-value', help='value for the attribute or range in the form "x y"'
                                                   ' (see --mjd)', nargs='*')
    query.add_argument('-d', '--dyn', action="store_true", help='Plot the dynamic spectrum')
    query.add_argument('-s', '--sec', action="store_true", help='Plot the secondary spectrum')
    # query.add_argument('-l', '--list', action="store_true", help='Only print a list of matching rows')
    query.add_argument('--csv', help='Print a csv list')

    parser.add_argument('database', help='filename of the database for reading and/or writing')
    parser.add_argument('-v', '--verbose', action="store_true", help='enable verbose mode')
    parser.add_argument('--debug', action="store_true", help='enable debug mode')
    args = parser.parse_args()
    if args.debug:
            print(args)
    return args


def main(args):
    """
    The main controller.
    :param args: arguments provided by the commandline
    """
    db = sqlite.DB(args.database, args.verbose, args.debug)
    if args.subcmd == 'ingest':
        db.ingest_data(args.files)
    elif args.subcmd == 'sql':
        for row in db.sql(args.sql):
            print(tuple(row))
    elif args.subcmd == 'query':
        attr_dict = {}
        if args.filename:
            attr_dict["filename"] = args.filename
        if args.pulsar:
            attr_dict["SOURCE"] = args.pulsar
        if args.mjd:
            attr_dict["MJD"] = args.mjd
        if args.freq:
            attr_dict["FREQ"] = args.freq
        rows = db.extract(attr_dict)

        # for row in rows:
        #     print(tuple(rows))
        print("Found {0} matching rows\n".format(len(rows)))
        # keys_dict = {'id': 5, 'filename': 30, 'ORIGIN': 10, 'MJD': 18, 'FREQ': 11, 'BW': 10}
        keys = ['id', 'filename', 'ORIGIN', 'MJD', 'FREQ', 'BW']
        widths = [4, 30, 23, 18, 12, 10]  # widths for the columns
        if len(rows) > 0:
            text, text2 = '', ''
            for key, value in zip(keys, widths):
                text += ' {0:{1}} |'.format(key, value)
                text2 += '{0}|'.format('-'*(value+2))  # makes the hline under the header +2 for the two spaces
            print(text[:-2])  # remove the last "\t|"
            # print('-'*(sum(widths)+len(widths)*3))
            print(text2[:-1])  # remove the last "|"
        csv = ''
        for key in rows[0][1].keys():
            if key != 'DATA':
                csv += '{0};'.format(key)
        csv += '\n'
        for row in rows:
            hdulist = row[0]
            db_header = row[1]
            if args.csv:
                for key in db_header.keys():
                    if key != 'DATA':
                        csv += '{0};'.format(db_header[key])
                csv += '\n'

            text = ''
            for key, value in zip(keys, widths):
                text += ' {0:{1}} |'.format(db_header[key], value)
            print(text[:-2])  # remove the last " |"
            # text = ''
            # for key in keys:
            #     text += '{0} = {1}, '.format(key, db_header[key])
            # print(text[:-2], '\n')
            if args.dyn:
                if not args.sec:  # only plot dynamic
                    dyn = computing.Dynamic(hdulist, db_header)
                    plotting.show_dyn(dyn)
            if args.sec:  # plot dynamic and secondary
                sec = computing.Secondary(hdulist, db_header)
                if args.dyn:
                    # plotting.show_image(sec.dyn)
                    plotting.show_dyn(sec)
                plotting.show_sec(sec)
            plotting.show()
        if args.csv:
            with open(args.csv, 'w') as f:
                f.write(csv)


if __name__ == '__main__':
    # for testing the computing and plotting
    # db = sqlite.DB('test.sqlite')
    # attrs = {'MJD': "0 60000"}
    # print(db.extract_data(attrs))
    args = parse_args()
    main(args)



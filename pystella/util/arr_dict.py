import csv

__author__ = 'bakl'


def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    :param dict_args:
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


def dict_save(dictionary, fname):
    """
    Save dict to file. Keys are column's names, values are column data
    :param dictionary:
    :type fname: the name of the file
    """
    with open(fname, 'wb') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(dictionary.keys())
        for row in zip(*dictionary.values()):
            # writer.writerow(list(row))
            writer.writerow(['{:12f}'.format(x) for x in row])
            # writer.writerow(['{:3.4e}'.format(x) for x in row])

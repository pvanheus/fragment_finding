#!/usr/bin/env python

import click
import sys
from bx.intervals.intersection import IntervalTree


class Interval(object):
    contig_id = start = end = strand = feature_id = None

    def __repr__(self):
        return 'Interval<contig={}, {}, {}, {}, ID={}>'.format(self.contig_id,
                                                               self.start, self.end, self.strand,
                                                               self.feature_id)


def parse_gff_attributes(attr_string):
    attributes = dict()
    attrs = attr_string.split(';')
    for attr in attrs:
        attr = attr.strip()
        if attr == '':
            continue
        (key, value) = attr.split('=')
        attributes[key] = value
    return attributes


def parse_gff_interval(infile, feature_type=None):
    count = 0
    for line in infile:
        count += 1
        if line.startswith('#'):
            continue
        fields = line.split('\t')
        assert len(fields) == 9, "Invalid GFF3 format at line {} in {}:\n{}".format(count, infile.name, line)
        # filter by feature type
        if feature_type is not None and fields[2] != feature_type:
            continue
        contig_id = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        attribute_string = fields[8]
        attributes = parse_gff_attributes(attribute_string)
        assert 'ID' in attributes, "Invalid GFF3 format, missing ID, on line {}:\n{}".format(count, line)
        interval = Interval()
        interval.contig_id = contig_id
        interval.start = start
        interval.end = end
        interval.strand = strand
        interval.feature_id = attributes['ID']
        yield interval

@click.command()
@click.argument('gene_gff3', type=click.File())
@click.argument('protein_gff3', type=click.File())
@click.argument('output_file', type=click.File(mode='w'), required=False, default=sys.stdout)
def print_fragmented_genes(gene_gff3, protein_gff3, output_file):
    # use of interval trees taken from this example:
    # http://informatics.malariagen.net/2011/07/07/using-interval-trees-to-query-genome-annotations-by-position/

    gene_trees = dict()
    genes = dict()
    for interval in parse_gff_interval(gene_gff3, feature_type='gene'):
        genes[interval.feature_id] = interval
        contig_tree = gene_trees.setdefault(interval.contig_id, dict())
        strand_tree = contig_tree.setdefault(interval.strand, IntervalTree())
        strand_tree.add(interval.start, interval.end, interval)
        gene_trees[interval.contig_id][interval.strand] = strand_tree

    fragment_set = set()
    for interval in parse_gff_interval(protein_gff3, feature_type='transcript'):
        genes = gene_trees[interval.contig_id][interval.strand].find(interval.start, interval.end)
        if len(genes) > 1:
            fragment_set.add(tuple([gene.feature_id for gene in genes]))

    count = 0
    for fragment_tuple in fragment_set:
        click.echo(fragment_tuple, file=output_file)
        count += len(fragment_tuple)
    click.echo(count, file=output_file)

if __name__ == '__main__':
    print_fragmented_genes()

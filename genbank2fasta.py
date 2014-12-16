#! /usr/bin/env python3
import argparse
from parsers import parse_gb, json_to_fasta, emit_fasta


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('gb_file', help="the GenBank file to convert")
  parser.add_argument('--tipdate',
                      action='store_true',
                      default=False,
                      help="Disallow sequences without collection dates")
  parser.add_argument('--filter',
                      nargs=argparse.REMAINDER,
                      default=[],
                      dest='filter_out',
                      help="keywords to disqualify sequences")
  args = parser.parse_args()
  gb_fh = open(args.gb_file, 'r')
  seqs_json = parse_gb(gb_fh.readlines())
  seqs_fasta = json_to_fasta(seqs_json, collection_date=args.tipdate,
      filter_out=args.filter_out)
  emit_fasta(seqs_fasta)

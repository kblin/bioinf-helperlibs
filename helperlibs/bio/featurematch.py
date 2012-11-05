#!/usr/bin/env python
#
# Find and make useable sequence features by locus tag.
# Copyright (C) 2008-2012 Kai Blin <kai.blin@biotech.uni-tuebingen.de>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

class FeatureMatch:
    """mRNA of a feature"""
    def __init__(self, feature, seq, direction, ul):
        self.feature = feature
        if direction == -1:
            self.direction = "reverse"
            self.dna = seq[ul:-ul].reverse_complement()
            self.long_dna = seq.reverse_complement()
            self.promotor_region = seq[-ul:].reverse_complement()
            self.mrna = seq.transcribe().reverse_complement()
        else:
            self.direction = "forward"
            self.dna = seq[ul:-ul]
            self.long_dna = seq
            self.promotor_region = seq[:ul]
            self.mrna = seq.transcribe()

        self.aas = self.dna.translate(to_stop=True)

    def get_fasta_header(self):
        ret = ">"
        if self.feature.qualifiers.has_key('locus_tag'):
            ret += self.feature.qualifiers['locus_tag'][0]
        elif self.feature.qualifiers.has_key('gene'):
            ret += self.feature.qualifiers['gene'][0]
        else:
            ret += "untagged"
        if self.feature.qualifiers.has_key('product'):
            ret += "|%s" % self.feature.qualifiers['product'][0]
        if self.feature.qualifiers.has_key('protein_id'):
            ret += "|%s" % self.feature.qualifiers['protein_id'][0]
        return ret

    def __str__(self):
        ret = "Feature:"
        if self.feature.qualifiers.has_key('locus_tag'):
            ret += "\n\tTag: %s" % self.feature.qualifiers['locus_tag'][0]
        elif self.feature.qualifiers.has_key('gene'):
            ret += "\n\tTag: %s" % self.feature.qualifiers['gene'][0]
        else:
            ret += "\n\tuntagged"
        ret += "\n\tStrand: %s" % self.direction
        if self.feature.qualifiers.has_key('product'):
            ret += "\n\tProduct: %s" % self.feature.qualifiers['product'][0]
        if self.feature.qualifiers.has_key('protein_id'):
            ret += "\n\tProtein ID: %s" % \
                self.feature.qualifiers['protein_id'][0]
        ret += "\n\tDNA: %s" % self.dna
        ret += "\n\tmRNA: %s" % self.mrna
        ret += "\n\tProtein: %s" % self.aas
        return ret

    def dna_fasta(self):
        """get feature DNA sequence in FASTA format"""
        ret = self.get_fasta_header()
        ret += "\n%s" % self.dna
        return ret

    def long_dna_fasta(self):
        """get feature DNA sequence including UTR in FASTA format"""
        ret = self.get_fasta_header()
        ret += "\n%s" % self.long_dna
        return ret

    def mrna_fasta(self):
        """get feature mRNA sequence in FASTA format"""
        ret = self.get_fasta_header()
        ret += "\n%s" % self.mrna
        return ret

    def protein_fasta(self):
        """get feature protein sequence in FASTA format"""
        ret = self.get_fasta_header()
        ret += "\n%s" % self.aas
        return ret

    def promotor_fasta(self):
        """get feature upstream UTR DNA sequence in FASTA format"""
        ret = self.get_fasta_header()
        ret += "\n%s" % self.promotor_region
        return ret

def find_features(seqs, locus_tag="all", utr_len=200):
    """Find features in sequences by locus tag"""
    found_features = []

    for seq_i in seqs:
        for feature in seq_i.features:
            if feature.type == "CDS" and (locus_tag == "all" or \
                    (feature.qualifiers.has_key('locus_tag') and \
                    feature.qualifiers['locus_tag'][0] == locus_tag)):
                start = max(0, feature.location.nofuzzy_start - utr_len)
                stop  = max(0, feature.location.nofuzzy_end + utr_len)
                feature_seq = seq_i.seq[start:stop]
                f_match = FeatureMatch(feature, feature_seq, feature.strand,
                                       utr_len)
                found_features.append(f_match)

    return found_features


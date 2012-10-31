#!/usr/bin/env python
#
# Thin convenience wrapper for BioPython's SeqIO code
# Copyright (C) 2012 Kai Blin <kai.blin@biotech.uni-tuebingen.de>
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

from Bio import SeqIO
from os import path

def parse(handle):
    dummy, ext = path.splitext(handle.name.lower())
    if ext in (".gbk", ".gb", ".genbank"):
        seq_iterator = SeqIO.parse(handle, "genbank")
    elif ext in (".embl"):
        seq_iterator = SeqIO.parse(handle, "embl")
    elif ext in (".fa", ".fasta", ".fna", ".faa"):
        seq_iterator = SeqIO.parse(handle, "fasta")
    else:
        handle.close()
        raise ValueError("Unknown file format '%s'." % ext)

    return seq_iterator

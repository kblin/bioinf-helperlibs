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
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO

def _get_seqtype_from_ext(handle):
    if isinstance(handle, basestring):
        name = handle
    elif hasattr(handle, 'name'):
        name = handle.name
    else:
        raise ValueError("Unknown datatype for handle!")
    dummy, ext = path.splitext(name.lower())
    if ext in (".gbk", ".gb", ".genbank", ".gbff"):
        return "genbank"
    elif ext in (".embl", ".emb"):
        return "embl"
    elif ext in (".fa", ".fasta", ".fna", ".faa", ".fas"):
        return "fasta"
    else:
        raise ValueError("Unknown file format '%s'." % ext)


def _guess_seqtype_from_file(handle):
    "Guess the sequence type from the file's contents"
    if isinstance(handle, basestring):
        handle = StringIO(handle)

    for line in handle:
        if not line.strip():
            continue
        if line.startswith('LOCUS'):
            return 'genbank'
        if line.startswith('ID'):
            return 'embl'
        if line.startswith('>'):
            return 'fasta'

    raise ValueError("Failed to guess format for input")

def parse(handle):
    seqtype = _get_seqtype_from_ext(handle)
    return SeqIO.parse(handle, seqtype)


def read(handle):
    seqtype = _get_seqtype_from_ext(handle)
    return SeqIO.read(handle, seqtype)


def write(*args, **kwargs):
    SeqIO.write(*args, **kwargs)

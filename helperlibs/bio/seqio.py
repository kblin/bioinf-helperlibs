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
    elif hasattr(handle, 'filename'):
        name = handle.filename
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
        if line.lstrip().split()[0] in ('LOCUS', 'FEATURES', 'source', 'CDS',
                                        'gene'):
            return 'genbank'
        if len(line) > 2 and line[:3] in ('ID ', 'FT '):
            return 'embl'
        if line.startswith('>'):
            return 'fasta'
    handle.seek(0)
    import string
    from Bio.Data import IUPACData as iupac
    all_input_letters = set(handle.read().lower())
    all_valid = set(string.digits)
    all_valid.update(set(iupac.protein_letters.lower()))
    all_valid.update(set(iupac.unambiguous_dna_letters.lower()))
    all_valid.update(set('- \n'))
    if all_valid.issuperset(all_input_letters):
        return 'fasta'

    raise ValueError("Failed to guess format for input")


def sanity_check_insdcio(handle, id_marker, fake_id_line):
    """Sanity check for insdcio style files"""
    found_id = False
    found_end_marker = False
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith(id_marker):
            found_id = True
            break
        if line.startswith('//'):
            found_end_marker = True
            break

    handle.seek(0)
    # We found an ID, file looks good.
    if found_id:
        return handle

    # If there's no ID and no end marker, just give up.
    if not found_end_marker:
        return handle

    # If we found an end marker but no ID, fake one.
    new_handle = StringIO()
    new_handle.write("%s\n" % fake_id_line)
    new_handle.write(handle.read())
    new_handle.seek(0)
    return new_handle


def sanity_check_embl(handle):
    """Sanity check EMBL format files"""
    id_marker = 'ID '
    fake_id_line = 'ID   DUMMY; SV 1; linear; DNA; STD; BCT; 1 BP.'
    return sanity_check_insdcio(handle, id_marker, fake_id_line)


def sanity_check_genbank(handle):
    """Sanity check GenBank format files"""
    id_marker = 'LOCUS '
    fake_id_line = 'LOCUS       DUMMY                      1 bp    DNA     linear   BCT 01-JAN-1970'
    return sanity_check_insdcio(handle, id_marker, fake_id_line)


def sanity_check_fasta(handle):
    """Sanity check FASTA files."""
    header_found = False
    for line in handle:
        if line.startswith('>'):
            header_found = True
            break

    handle.seek(0)

    if header_found:
        return handle

    fake_header_line = ">DUMMY"
    new_handle = StringIO()
    new_handle.write("%s\n" % fake_header_line)
    new_handle.write(handle.read())
    new_handle.seek(0)
    return new_handle


def parse(handle, robust=False):
    seqtype = _get_seqtype_from_ext(handle)
    if robust:
        if seqtype == "embl":
            handle = sanity_check_embl(handle)
        elif seqtype == "genbank":
            handle = sanity_check_genbank(handle)
        elif seqtype == "fasta":
            handle = sanity_check_fasta(handle)

    return SeqIO.parse(handle, seqtype)


def read(handle, robust=False):
    seqtype = _get_seqtype_from_ext(handle)
    if robust:
        if seqtype == "embl":
            handle = sanity_check_embl(handle)
        elif seqtype == "genbank":
            handle = sanity_check_genbank(handle)
        elif seqtype == "fasta":
            handle = sanity_check_fasta(handle)
    return SeqIO.read(handle, seqtype)


def write(*args, **kwargs):
    SeqIO.write(*args, **kwargs)

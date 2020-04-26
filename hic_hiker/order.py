#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
fixing of ordering
"""
import numpy as np
import pandas as pd

from .contigs import Contigs
from .load import get_contacts_mnd
from .load_3ddna import Assembly
from .prob import infer_from_longest_contig2

def load_all(fasta_filename, mnd_filename, asm_filename):
    # load all files
    contigs = Contigs(fasta_filename=fasta_filename)
    df = get_contacts_mnd(
            contigs=contigs,
            mnd_filename=mnd_filename
            )
    asm = Assembly(asm_filename=asm_filename, contigs=contigs)
    layout = asm.get_layout()
    contigs2, df2 = asm.update_all(contigs, df)
    estimator, raw_estimator = infer_from_longest_contig2(df2, contigs2, K=75000, K0=10**3.5)

    # convert them into pandas.DataFrame


    return df_contigs, df_contacts, df_layout, estimator


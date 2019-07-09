#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
import unittest
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
import contigs

class TestContigs(unittest.TestCase):
    def setUp(self):
        c = contigs.Contigs('../small_mock/contigs.fasta')
        print(c.names)
        self.c = c
        self.x = 10
    def test_hoge(self):
        x = self.x
        self.assertEqual(x, 10)

if __name__ == '__main__':
    unittest.main()

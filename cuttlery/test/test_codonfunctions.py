#!/usr/bin/env python

# cuttlery - 'Codon Usage Table Tools-lery.'
# Copyright (c) 2016-2017 Darrin T. Schultz. All rights reserved.
#
# This file is part of cuttlery.
#
# cuttlery is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cuttlery is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with cuttlery.  If not, see <http://www.gnu.org/licenses/>.

"""@author Darrin Schultz
This class tests the classes and methods for codonfunctions.py
"""

import unittest
import cuttlery.codonfunctions as ccf
import re
import os
import logging

class codonfunctions_test_case(unittest.TestCase):
    """Tests that the functions in codonfunctions work properly
    """
    def setUp(self):
        self.pwd = os.path.dirname(os.path.abspath(__file__))
        self.ideal_results = {'adh': '/Users/darrin/git/cuttlery/cuttlery/test/fasta_test_dir/adh.fa',
                         'COI': '/Users/darrin/git/cuttlery/cuttlery/test/fasta_test_dir/subdir1/COI.fasta'}

    def test_fail_with_int(self):
        #should fail with an int
        self.assertRaises(Exception,
                          ccf.fasta_dir_to_gene_filelist,
                          1)

    def test_fail_with_float(self):
        #should fail with an float
        self.assertRaises(
            Exception,
            ccf.fasta_dir_to_gene_filelist,
            2.1)

    def test_fail_with_boolean(self):
        #should fail with a boolean
        self.assertRaises(
            Exception,
            ccf.fasta_dir_to_gene_filelist,
            True)
    def test_fail_with_empty_list(self):
        # should fail with an empty list
        self.assertRaises(
            IOError,
            ccf.fasta_dir_to_gene_filelist,
            [])

    def test_work_with_single_directory_path(self):
        """Make sure that the program returns all directories properly
        when searching recursively through files
        """
        filelist = ccf.fasta_dir_to_gene_filelist(self.pwd)
        self.assertEqual(filelist, self.ideal_results)

    def test_work_with_multiple_overlapping_directories(self):
        """Make sure the method works when it is fed multple directories
        that are overlapping
        """
        pathlist = [self.pwd,
                    os.path.join(self.pwd,"fasta_test_dir"),
                    ps.path.join(self.pwd,"fasta_test_dir/subdir1")] 
        filelist = ccf.fasta_dir_to_gene_filelist(pathlist)
        self.assertEqual(filelist, self.ideal_results)

if __name__ == '__main__':
    unittest.main()

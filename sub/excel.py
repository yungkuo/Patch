# -*- coding: utf-8 -*-
"""
Created on Fri Oct 24 21:40:42 2014

@author: Philip
"""

import numpy as np
import xlrd


def openexcel(filename):
    workbook = xlrd.open_workbook(filename)
    worksheet = workbook.sheet_by_name('Sheet1')
    #num_rows = worksheet.nrows - 1
    #num_cells = worksheet.ncols - 1
    #curr_row = -1

    cell_value=np.zeros((worksheet.nrows, worksheet.ncols))
    for col in range(worksheet.ncols):
        for row in range(worksheet.nrows):
            cell_value[row, col] = worksheet.cell_value(row, col)

    return cell_value

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib as plt
from scipy import ndimage as ndi
from skimage import data, io, filters, measure, feature, exposure
from skimage.measure import regionprops
from skimage.morphology import watershed

# ----------------------------------------------------------

PRODUCE_IMAGES = False
# ----------------------------------------------------------

def get_nucleus_gfp(im_gfp, im_dna, border_grayscale_threshold):
    # custome filter: cut out only the "whiter" stain-affected areas
    # then paint the nucleus area with gray-scale values from the GFP image (channel1)
    imf = np.where(im_dna > border_grayscale_threshold, im_gfp, 1)    
    return imf

def get_cytoplasm_gfp(im_gfp, im_dna, border_grayscale_threshold):
    imf = np.where(im_dna <= border_grayscale_threshold, im_gfp, 1)    
    return imf
    
def get_image_gfp_points(im, border_grayscale_threshold):
    imf = im[im >= border_grayscale_threshold]
    return imf

def segment_well(image, intensity_image):
    # segment -----------------------------
    edges = filters.sobel(image)
    
    markers = np.zeros_like(image)
    markers[image < 5] = 1
    markers[image > 25] = 2
    
    # segment
    segs = watershed(edges, markers) #gets segments valued by : {0,1,2}
    segs_binary = ndi.binary_fill_holes(segs - 1)
    
    # count the segments
    labeled_segs, n_labels = ndi.label(segs_binary)    
    props = regionprops(labeled_segs, intensity_image)
    
    return props

def get_well_cyto_gfps(props, border_grayscale_threshold):
    min_box_size = (15, 15) # ***
    well_gfps = []
    cell_labels = []
    for prop in props:
        if prop.image.shape > min_box_size: #ignore incorrent very small segments
            compartment_img = prop.intensity_image
            fgp_points = get_image_gfp_points(compartment_img, border_grayscale_threshold)
            if len(fgp_points) > 0:
                gfp = np.mean(fgp_points)
                well_gfps.append(gfp)
                cell_labels.append(prop.label)
    return well_gfps, cell_labels

def get_well_nuc_gfps(props, border_grayscale_threshold, valid_cell_lables):
    well_gfps = []
    for prop in props:
        if prop.label in valid_cell_lables:
            compartment_img = prop.intensity_image
            fgp_points = get_image_gfp_points(compartment_img, border_grayscale_threshold)
            if len(fgp_points) > 0:
                gfp = np.mean(fgp_points)
                well_gfps.append(gfp)
            else:
                well_gfps.append(0.0)
    return well_gfps
            
# ----------------------------------------------
# ----------------------------------------------
    
gfp_nucleus_list = []
gfp_cytoplasm_list = []
gfp_total_list = []
dose_list = []
row_list = []
col_list = []
cell_list = []

plate_doses = np.loadtxt("./doses_BBBC013_v1_platemap_all.txt")

rows = ['A', 'B', 'C', 'D', 'E', 'F' ,'G', 'H']
row_index = 0
col = 1
plat_no = 1
for row_index in range(0,8):    #0,8 (A to H)
    for col in range(1,13):     #1,13
        row = rows[row_index]
        plate_name = "{:02}-{:}-{:02}.BMP".format(plat_no, row, col)
        filepath1 = './BBBC013_v1_images_bmp/' + "Channel1-" + plate_name
        filepath2 = './BBBC013_v1_images_bmp/' + "Channel2-" + plate_name
        print (plate_name)
        im_gfp = io.imread(filepath1) 
        im_dna = io.imread(filepath2)
        
        # Find GFP of each compartment (nucleus, cytoplasm) --------
        # GFP intensty values based on the grayscale intensities of pixels
        nucleus_border_grayscale_threshold = 25 # 20..30 seems to be the sweet spot by looking at the grayscale histogram
        im_cyto_gfp = get_cytoplasm_gfp(im_gfp, im_dna, nucleus_border_grayscale_threshold)
        im_nucleus_gfp = get_nucleus_gfp(im_gfp, im_dna, nucleus_border_grayscale_threshold)
        
        if PRODUCE_IMAGES:
            image_filename1 = './output_pics/{:}-nucleus-with-gfp.png'.format(plate_name)
            image_filename2 = './output_pics/{:}-cytoplasm-with-gfp.png'.format(plate_name)
            io.imsave(image_filename1, im_nucleus_gfp)
            io.imsave(image_filename2, im_cyto_gfp)

        # get other values -------
        dose = plate_doses[plat_no-1]
                   
        # Segment the well -----------------------
        border_grayscale_threshold = nucleus_border_grayscale_threshold
        # CYTOPLASMs:
        cyto_props = segment_well(im_cyto_gfp, im_cyto_gfp)
        cyto_gfps, valid_cell_labels = get_well_cyto_gfps(cyto_props, border_grayscale_threshold)
        # Count the cells
        for cell_no in range(1,len(cyto_gfps)+1):
            cell_list.append(cell_no)
            row_list.append(row)
            col_list.append(col)
            dose_list.append(dose)
        # Cytos
        for segment_gfp in cyto_gfps:            
            gfp_cytoplasm_list.append(segment_gfp)                        
                    
        # Nuclei
        nuc_props = segment_well(im_cyto_gfp, im_nucleus_gfp)
        nuc_gfps = get_well_nuc_gfps(nuc_props, border_grayscale_threshold, valid_cell_labels)
        for segment_gfp in nuc_gfps:            
            gfp_nucleus_list.append(segment_gfp)                        
        
        plat_no += 1

# Calculate Ratio of GFP in Nucleus VS. GFP in Cytoplasm
gfp_ratio_list = np.array(gfp_nucleus_list)/np.array(gfp_cytoplasm_list)

# Save values to CSV
pd.DataFrame(data={'row':row_list, 'column':col_list, 'Dose':dose_list, 'cell_no':cell_list, 'GFP_in_cytoplasm':gfp_cytoplasm_list, 'GFP_in_nucleus':gfp_nucleus_list, 'GFP_nuc_to_cyto_ratio':gfp_ratio_list}).to_csv("gfp_cells_complete.csv")

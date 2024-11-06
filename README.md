# Development-of-a-tool-for-identifying-and-selecting-rice-grains

## Project Overview

This repository provides code for analyzing binary images of rice grains. The code uses **DGtal** library functionalities to process these images, extract digital objects representing rice grains, and calculate geometric properties, such as area, perimeter, and circularity. It also includes steps to visualize and save the output.

## Dependencies

The code requires the following:
- **DGtal**: Digital geometry library for topological and geometric operations.
- **LibBoard**: Library for creating and manipulating vector images in EPS format.

## Code Description

### 1. Import Libraries and Define Types
The code imports essential **DGtal** headers for image processing, digital topology, and geometry handling.

### 2. Main Workflow
The `main` function processes a list of rice grain images and performs the following steps:

1. **Image Loading and Digital Set Creation**  
   Reads binary PGM images, converts them into a digital set, and visualizes it.

2. **Connected Component Extraction**  
   Identifies and separates individual rice grains as connected components.

3. **Boundary Extraction and Visualization**  
   Extracts and saves the boundaries of each rice grain.

4. **Polygonization**  
   Converts boundaries to polygonal representations using DSS and saves the results.

5. **Area Calculation**  
   Calculates areas of each rice grain using two methods:
   - Counting the number of pixels in the digital set (2-cells).
   - Shoelace formula on the polygon vertices.

6. **Perimeter Calculation**  
   Calculates perimeters using:
   - 1-cell boundary pixels.
   - Polygon segment lengths.

7. **Circularity Calculation**  
   Computes circularity for each rice grain as a measure of roundness.

## Usage

1. **Compile the Code**  
   Use a C++ compiler compatible with DGtal and LibBoard.

2. **Run the Code**  
   Run the compiled program to process images from the `./RiceGrains/` folder.


### Output
The program outputs EPS files for each step and prints properties (area, perimeter, circularity) to the console.


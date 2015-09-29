#ifndef __IMAGE_H
#define __IMAGE_H

#include <cstdio>
#include <ctime>
#include <cstring>
#include <cmath>
#include <assert.h>
#include "graph.h"
#include "kovtun.h"
#include "aexpand.h"
#include "TRWBP.h"

/**** Assumes a grid (image-like) graph ****/
class image
{
	public:
	int width;           // Grid width
	int height;          // Grid height
	int num_nodes;       // Number of nodes in the grid
	char filename[100];  // Image filename (assumes a PPM file)
	int num_labels;      // Number of labels

	Energy *eng;               // Stores the unary and pairwise costs
	unsigned char **label_map; // Label assignments to pixels
	unsigned char **buf;       // Pixel colour values

	image(const char *file, int nlabel);
	~image();

	/**** multi-label solver functions ****/
	void kovtun(bool do_aexpand=false);   // Computes the partially optimal solution using Kovtun's method.
	                                      // The remaining nodes can be solved using alpha expansion if do_expand = true.
	void trw(bool use_projection=false);  // use_projection = true : Partially optimal solution computation + TRW/BP on remaining nodes
	                                      // use_projection = false: Standard TRW/BP
	/*************************************/

	/** Image file read/write functions **/
	void readPPM();
	void readPGM();
	void readPPMWO();
	void writePPM(const char *file);
	void writePGM(const char *file);
};

#endif

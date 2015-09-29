#include "image.h"

#define OBJ 0
#define BKG 1
#define ACTIVE true
#define INACTIVE false

/************************************************************************/
/***************** See image.h for function description *****************/
/************************************************************************/

image::image(const char *file, int nlabel)
{
	int i, j;

	strcpy(filename, file);
	readPPM();
	num_labels = nlabel;
	num_nodes = width*height;

	label_map = new unsigned char* [height];
	for(i=0; i<height; i++)
		label_map[i] = new unsigned char [width];

	eng = new Energy(num_labels);
	eng->nvar = num_nodes;
	eng->npair= 4*width*height - 3*(width+height) + 2; // Using an 8-neighbourhood
	eng->allocateMemory();

	/**
	 ** Assign the unary (eng->unaryCost) and pairwise costs (eng->pairCost).
	 ** Note that pairwise costs are indexed by eng->pairIndex, which is also to be initialized.
	 **/
}

image::~image()
{
	int i;

	for(i=0; i<height; i++)
	{
		delete buf[i];
		delete label_map[i];
	}
	delete buf;
	delete label_map;
	delete eng;
}

void image::kovtun(bool do_aexpand)
{
	int i, j;
	unsigned char lbl;
	double time_np;
	Grapht** graph;

	for(i=0; i<height; i++)
	for(j=0; j<width; j++)
		label_map[i][j] = NOLABEL;

	graph = new Grapht* [num_labels];
	for(i=0; i<num_labels; i++) graph[i] = new Grapht(num_nodes, eng->npair);
	Kovtun* solver = new Kovtun(eng, graph, 1); // 1 => use multiple graphs and the Reuse method (see Alahari et al., 2008)

	/* Projecting the energy function (Reduce method, see Alahari et al., 2008) */
	for(i=0; i<num_labels; i++)
	{
		time_np = solver->findPersistent(i);
		eng->Project(solver->multiSolution);
	}

	if(do_aexpand)
	{
	/* Solve the remaining nodes using the alpha expansion method */
		Aexpand *solvera;
		eng->Project(solver->multiSolution);
		solvera = new Aexpand(eng, graph);
		time_np = solvera->minimize(solver->nodes);

		int count=0;
		for(i=0; i<height; i++)
		for(j=0; j<width; j++)
		{
			lbl = solver->multiSolution[i*width+j];
			if(lbl != NOLABEL)
				label_map[i][j] = lbl;
			else
			{
				label_map[i][j] = solvera->label_map[eng->activeU[count]];
				count++;
			}
		}
		delete solvera;
	}
	else
	{
	/* Or compute only the partially optimal solution */
		for(i=0; i<height; i++)
		for(j=0; j<width; j++)
		{
			lbl = solver->multiSolution[i*width+j];
			if(lbl != NOLABEL)
				label_map[i][j] = lbl;
		}
	}

	for(i=0; i<num_labels; i++)
		delete graph[i];
	delete [] graph;
	delete solver;
}

void image::trw(bool use_projection)
{
	int i, j;
	unsigned char lbl;
	Energy* n_eng;

	for(i=0; i<height; i++)
	for(j=0; j<width; j++)
		label_map[i][j] = NOLABEL;

	if(use_projection)
	{
	/* Compute the partially optimal solution and solve the remaining nodes using TRW and BP */
		Grapht* graph  = new Grapht(num_nodes, eng->npair);
		Kovtun* solver = new Kovtun(eng, graph, 0); // 0 => use a single graph
		double time_np;

		for(i=0; i<num_labels; i++)
		{
			time_np = solver->findPersistent(i);
			eng->Project(solver->multiSolution);
			graph->reset();
		}
		n_eng = eng->Projection(solver->multiSolution);

		TRWBP *solvert = new TRWBP(n_eng);
		solvert->minimize(true, 70, 70); // 70 iterations each of BP and TRW

		int count=0;
		for(i=0; i<height; i++)
		for(j=0; j<width; j++)
		{
			lbl = solver->multiSolution[i*width+j];
			if(lbl != NOLABEL)
				label_map[i][j] = lbl;
			else
			{
				// solvert->solutionT contains TRW result and 
				// solvert->solutionB contains BP result
				label_map[i][j] = solvert->solutionT[count];
				solver->multiSolution[i*width+j] = label_map[i][j]; // To compute energy
				count++;
			}
		}
		delete graph;
		delete solver;
		delete solvert;
		delete n_eng;
	}
	else
	{
	/* Or solve using standard TRW and BP only */
		n_eng = eng;
		TRWBP *solver = new TRWBP(n_eng);
		solver->minimize(true, 70, 70); // 70 iterations each of BP and TRW

		// solvert->solutionT contains TRW result and 
		// solvert->solutionB contains BP result
		for(i=0; i<height; i++)
		for(j=0; j<width; j++)
			label_map[i][j] = solver->solutionT[i*width+j];
		delete solver;
	}
}


/* Functions to read/write images */
void image::readPPM()
{
	int i,j;
	int Width,Height;

	FILE *fp = fopen(filename, "r");

	fscanf(fp,"P6\n");
	while(1)
	{
		char s[80];
		fgets(s,80,fp);
		if(s[0] >= '0' && s[0] <= '9')
		{
			sscanf(s,"%d %d", &width, &height);
			break;
		}
	}
	fscanf(fp, "255\n");
	
	Width = 3*(width);
	Height= height;

	buf = new unsigned char* [Height];
	for(i=0; i<Height; i++)
		buf[i] = new unsigned char [Width];

	for(i=0; i<Height; i++)
		for(j=0; j<Width; j++)
			fscanf(fp, "%c", &buf[i][j]);
	
	fclose(fp);
}

void image::readPGM()
{
	int i,j;
	int Width,Height;

	FILE *fp = fopen(filename, "r");

	fscanf(fp,"P5\n");
	while(1)
	{
		char s[80];
		fgets(s,80,fp);
		if(s[0] >= '0' && s[0] <= '9')
		{
			sscanf(s,"%d %d", &width, &height);
			break;
		}
	}
	fscanf(fp, "255\n");
	
	buf = new unsigned char* [height];
	for(i=0; i<height; i++)
		buf[i] = new unsigned char [width];

	for(i=0; i<height; i++)
		for(j=0; j<width; j++)
			fscanf(fp, "%c", &buf[i][j]);
	
	fclose(fp);
}

void image::readPPMWO()
{
	int i,j;
	int Width,Height;

	FILE *fp = fopen(filename, "r");

	fscanf(fp,"P6\n");
	while(1)
	{
		char s[80];
		fgets(s,80,fp);
		if(s[0] >= '0' && s[0] <= '9')
		{
			sscanf(s, "%d%d", &width, &height);
			break;
		}
	}
	fscanf(fp, "255\n");

	Width = 3*(width);
	Height= height;

	for(i=0; i<Height; i++)
		for(j=0; j<Width; j++)
			fscanf(fp, "%c", &buf[i][j]);

	fclose(fp);
}

void image::writePPM(const char *file)
{
	FILE *fp = fopen(file, "w");
	int i,j;

	fprintf(fp, "P6\n%d %d\n255\n", width, height);
	width *= 3;
	for(i=0; i<height; i++)
		for(j=0; j<width; j++)
			fprintf(fp, "%c", buf[i][j]);
	fclose(fp);
	width /= 3;
}

void image::writePGM(const char *file)
{
	FILE *fp = fopen(file, "w");
	int i,j;

	fprintf(fp, "P5\n%d %d\n255\n", width, height);
	for(i=0; i<height; i++)
		for(j=0; j<width; j++)
			fprintf(fp, "%c", buf[i][j]);
	fclose(fp);
}


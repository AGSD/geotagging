#include "include/pcp.h"
#include "include/spatial_helpers.h"
#include <fstream>
#include <cstdio>

#define DBG(x)

unsigned int numIter = 0;
inline bool validLocation(location &l){
	double lat = l.lat;
	double lng = l.lng;
	if((-90.0<=lat && lat<=90.0) && (-180.0<=lng && lng<=180.0))
		return true;
	return false;
}

struct GT{
	location **locationLists;
    unsigned int *validIndex;
    location *nodeLocation;
    bool *locKnown;
    
    GT(graph<location> &G, location *nodeLocationInit, bool *locKnownIn){	//constructor
		int n = G.numVertex;
		locKnown = locKnownIn;
    	validIndex = new unsigned int[n+1];
    	locationLists = new location*[n+1];
    	nodeLocation = nodeLocationInit;
    	for(unsigned int i=1; i<=n; ++i){
    		if(!locKnown[i])
	    		locationLists[i] = new location[G.inDeg[i]];
	    	else
	    		locationLists[i] = NULL;
     		validIndex[i] = 0;
     	}
     	DBG(printf("GT init completed\n");)
    }

    inline location scatterFunc (unsigned int node)
    {
        return nodeLocation[node];
    }

    inline bool reInit(unsigned int node)
    {
        return true;
    }

    inline bool gatherFunc (location updateVal, unsigned int destId)
    {	
    	if(!locKnown[destId] && validLocation(updateVal))
	        locationLists[destId][validIndex[destId]++] = updateVal;
	        
        return true;
    }  
    
    inline bool apply(unsigned int node)
    {
    	if(!locKnown[node])
	    	nodeLocation[node] = spatial_center(locationLists[node], validIndex[node]);
        return true;
    } 

};

int main(int argc, char** argv)
{
	if(argc<3){
		printf("Usage: ./geotagging.o <graphFilename> <locationFilename> -t <numThreads(optional)>  -rounds <#rounds(default 3)>\n");
		return 0;
	}
    graph<location> G;
    initialize(&G, argc, argv);
    initBin<location>(&G);
	unsigned int n = G.numVertex;
    DBG(printf("N = %u\n",n);)
    location *nodeLocation = new location[n+1];
    bool *locKnown = new bool[n+1]();
    
    DBG(printf("nodeLocation[1] = %lf %lf\n",nodeLocation[1].lat,nodeLocation[1].lng);)
    unsigned int initFrontierSize = 0;
    unsigned int *tmpFrontier = new unsigned int [n+1];
	            
    //reading in locations
    char* locationFile = argv[2];
    ifstream f;
    f.open(locationFile);
    unsigned int nodeId;
    double lat,lng;

    while(f>>nodeId>>lat>>lng){
    	nodeLocation[nodeId].lat = lat;
	   	nodeLocation[nodeId].lng = lng;
	   	tmpFrontier[initFrontierSize++] = nodeId;
	   	locKnown[nodeId] = true;
    	//dbg(printf("nodeId = %u, lat = %lf, lng = %lf\n",nodeId,lat,lng);)
    }
    f.close();
	dbg(printf("nodeId = 4, lat = %lf, lng = %lf, outDeg = %u\n",nodeLocation[4].lat,nodeLocation[4].lng,G.outDeg[4]);)
    unsigned int * initFrontier = new unsigned int [initFrontierSize];
    
    DBG(printf("Frontier size = %u",initFrontierSize);)
    for (int i=0; i<initFrontierSize; i++)
        initFrontier[i] = tmpFrontier[i];

    delete tmpFrontier;

    struct timespec start, end, half;
    float time;
    
    numIter = 0;
    int ctr =0;
	
	DBG(printf("entering loop\n");)
	
    while(ctr < G.rounds){
        loadFrontier(&G, initFrontier, initFrontierSize);

        if( clock_gettime(CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}
		DBG(printf("entering inner loop\n");)
        while(numIter<1)
        {
        	DBG(printf("inside inner loop\n");)
            pcpm<location>(&G, GT(G,nodeLocation,locKnown));
            numIter++;
        }

        if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");}
        time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;
        printf("geotagging, %d, %s, %lf\n", NUM_THREADS, argv[1], time);
        ctr++;
    }
    
    ofstream o;
    o.open("output.txt");
    
    for(int i=1; i<=n; ++i){
    	o<<i<<" "<<nodeLocation[i].lat<<" "<<nodeLocation[i].lng<<'\n';
    }
    
    o.close();
    
    printf("\n");

    return 0;
}




#include "../include/mesh.h"

Mesh* newMesh (int argc, char *argv[])
{
    printf("--------------------------------------------------------------------------------------------\n");
    printf("[!] Reading VTK file '%s' !\n",argv[1]);
    Mesh *mesh = (Mesh*)malloc(sizeof(Mesh));
    mesh->h = setTypeCell(argv[3]);
    Graph *g = readPurkinjeNetworkFromFile(argv[1]);
    GraphToMesh(mesh,g);
    writeLevelToFile(mesh,g);
    //printGraph(g);
    
    #ifdef DEBUG
    writeMeshInfo(mesh);
    #endif
    
    printf("--------------------------------------------------------------------------------------------\n");
    return mesh;
}

// Return the length of cell based on name of cell given by the input
double setTypeCell (char cName[])
{
    size_t found;
    string str(cName);
    found = str.find("alien");
    if (found != string::npos)
        return ALIEN_DX;
    found = str.find("dog");
    if (found != std::string::npos)
        return DOG_DX;
    found = str.find("orc");
    if (found != string::npos)
        return ORC_DX;
    found = str.find("pig");
    if (found != string::npos)
        return PIG_DX;

    found = str.find("test1");
    if (found != string::npos)
        return TEST1_DX;
    found = str.find("test2");
    if (found != string::npos)
        return TEST2_DX;
    found = str.find("test3");
    if (found != string::npos)
        return TEST3_DX;
    found = str.find("test4");
    if (found != string::npos)
        return TEST4_DX;
    found = str.find("test5");
    if (found != string::npos)
        return TEST5_DX;
    found = str.find("test6");
    if (found != string::npos)
        return TEST6_DX;
    return 0;
}

// Calculate the direction vector from node 1->2 and return its unitary vector
void calcDirection (double d_ori[], Node *n1, Node *n2)
{
    double norm;
    d_ori[0] = n2->x - n1->x; d_ori[1] = n2->y - n1->y; d_ori[2] = n2->z - n1->z;
    norm = sqrt(pow(d_ori[0],2)+pow(d_ori[1],2)+pow(d_ori[2],2));
    for (int i = 0; i < 3; i++) d_ori[i] /= norm;
}

// Map the ids from the graph to the ids of mesh 
void GraphToMesh (Mesh *mesh, Graph *g)
{
    int N = g->total_nodes;
    mesh->map_graph_elem = (int*)calloc(N,sizeof(int));
    // Construir o nodo 0 primeiro
    Node *ptr = g->listNodes;
    Point point; point.x = ptr->x; point.y = ptr->y; point.z = ptr->z; point.d = D;
    mesh->points.push_back(point);
    // Busca em profundidade
    DFS(mesh,ptr,0);
}

// Depth-First-Search
void DFS (Mesh *mesh, Node *u, int lvl)
{
    // Diameter decrease rule
    double d = D*pow(DELTA,lvl);
    //printf("Node %d || Level %d || Diameter %lf\n",u->id,lvl,d);
    Edge *v = u->edges;
    while (v != NULL)
    {
        growSegment(mesh,u,v,d);
        DFS(mesh,v->dest,lvl+1);
        v = v->next;
    }
    mesh->nPoints = mesh->points.size();
    mesh->nElem = mesh->elements.size();
}

void growSegment (Mesh *mesh, Node *u, Edge *v, double diam)
{
    double d_ori[3], d[3];
    double w = v->w;
    int nElem = w / mesh->h;
    int id_source = mesh->map_graph_elem[u->id];
    calcDirection(d_ori,u,v->dest);
    d[0] = u->x; d[1] = u->y; d[2] = u->z;
    printf("Node %d will grow %d points\n",u->id,nElem);
    // Grow the number of elements of size 'h' until reach the size of fiber
    for (int k = 1; k <= nElem; k++)
    {
        double x, y, z;
        x = d[0] + d_ori[0]*mesh->h*k; y = d[1] + d_ori[1]*mesh->h*k; z = d[2] + d_ori[2]*mesh->h*k;
        Point point; point.x = x; point.y = y; point.z = z; point.d = diam;
        mesh->points.push_back(point);
        Element elem; elem.left = id_source; elem.right = mesh->points.size()-1;
        mesh->elements.push_back(elem);
        id_source = mesh->points.size()-1;
    }
    // Save the of the last inserted node. In case this node generate offsprings
    mesh->map_graph_elem[v->id] = id_source;
}   

// Print information about the mesh to a file
void writeMeshInfo (Mesh *mesh)
{
    FILE *outFile = fopen("infoMesh.txt","w+");
    fprintf(outFile,"---------------------------- MESH INFO --------------------------------------------\n");
    fprintf(outFile,"Number of points = %d\n",mesh->nPoints);
    fprintf(outFile,"Number of elements = %d\n",mesh->nElem);
    fprintf(outFile,"h = %lf\n\n",mesh->h);
    fprintf(outFile,"-> POINTS\n");
    for (int i = 0; i < (int)mesh->points.size(); i++)
        fprintf(outFile,"Point %d = (%lf,%lf,%lf)\n",i,mesh->points[i].x,mesh->points[i].y,mesh->points[i].z);
    fprintf(outFile,"-> ELEMENTS\n");
    for (int i = 0; i < (int)mesh->elements.size(); i++)
        fprintf(outFile,"Element %d = (%d,%d)\n",i,mesh->elements[i].left,mesh->elements[i].right);
    fprintf(outFile,"-----------------------------------------------------------------------------------\n");
    fclose(outFile);
}

// Print information about the mesh to the standart output
void printMeshInfo (Mesh *mesh)
{
    printf("---------------------------- MESH INFO --------------------------------------------\n");
    printf("Number of points = %d\n",mesh->nPoints);
    printf("Number of elements = %d\n",mesh->nElem);
    printf("h = %lf\n\n",mesh->h);
    printf("-> POINTS\n");
    for (int i = 0; i < (int)mesh->points.size(); i++)
        printf("Point %d = (%lf,%lf,%lf)\n",i,mesh->points[i].x,mesh->points[i].y,mesh->points[i].z);
    printf("-> ELEMENTS\n");
    for (int i = 0; i < (int)mesh->elements.size(); i++)
        printf("Element %d = (%d,%d)\n",i,mesh->elements[i].left,mesh->elements[i].right);
    printf("-----------------------------------------------------------------------------------\n");
}

// Write the mesh to a file (.msh)
void writeMeshToFile (Mesh *mesh, char *filename)
{
    if (mesh != NULL)
    {
        //changeExtension(filename);
        printf("[!] Mesh file will be saved in :> %s\n",filename);
        FILE *file = fopen(filename,"w+");
        if (!file) printf("[-] ERROR! Opening file '%s'\n",filename);
        fprintf(file,"%d %d %lf\n",mesh->nElem,mesh->nPoints,mesh->h);
        for (int i = 0; i < mesh->nPoints; i++)
        {
            // Malha com a informacao adicional do diametro nos pontos
            #ifdef DIAMETER
            fprintf(file,"%lf %lf %lf %lf\n",mesh->points[i].x,mesh->points[i].y,mesh->points[i].z,mesh->points[i].d);
            // Malha padrao so com a informacao das coordenadas dos pontos
            #else
            fprintf(file,"%lf %lf %lf\n",mesh->points[i].x,mesh->points[i].y,mesh->points[i].z);
            #endif
        }
        for (int i = 0; i < mesh->nElem; i++)
            fprintf(file,"%d %d\n",mesh->elements[i].left,mesh->elements[i].right);
        fclose(file);
    }   
}

// Write the mesh to a file (.vtk)
void writeMeshToVTK (Mesh *mesh, const char *filename)
{
    if (mesh != NULL)
    {
        //changeExtension(filename);
        printf("[!] Mesh file will be saved in :> %s\n",filename);
        FILE *outFile = fopen(filename,"w+");
        fprintf(outFile,"# vtk DataFile Version 3.0\n");
        fprintf(outFile,"Purkinje\nASCII\nDATASET POLYDATA\n");
        fprintf(outFile,"POINTS %d float\n",(int)mesh->points.size());
        for (int i = 0; i < (int)mesh->points.size(); i++)
            fprintf(outFile,"%.8lf %.8lf %.8lf\n",mesh->points[i].x,mesh->points[i].y,mesh->points[i].z);
        fprintf(outFile,"LINES %d %d\n",(int)mesh->elements.size(),(int)mesh->elements.size()*3);
        for (int i = 0; i < (int)mesh->elements.size(); i++)
            fprintf(outFile,"2 %d %d\n",mesh->elements[i].left,mesh->elements[i].right);
        fclose(outFile);
    }   
}

// Change the extension of a file
void changeExtension (char *filename)
{
    int N = strlen(filename);
    filename[N-3] = 'm';
    filename[N-2] = 's';
    filename[N-1] = 'h';
    printf("[!] Mesh file will be saved in :> %s\n",filename);
}

// Template to the sort function of the 'map' structure
template <typename T1, typename T2>
struct less_second
{
    typedef pair<T1,T2> type;
    bool operator ()(type const &a, type const &b) const
    {
        return a.second < b.second;
    }
};

// Run a Breath-First-Seach and write the levels of the network
void writeLevelToFile (Mesh *mesh, Graph *g)
{
    FILE *lvlFile = fopen("level.txt","w+");
    map<int,int> dist = BFS(g,0);
    // Copy the map to a vector
    vector< pair<int,int> > levels(dist.begin(),dist.end());
    // Sort by level
    sort(levels.begin(),levels.end(),less_second<int,int>());
    // Print the level
    int maxLevel = levels[levels.size()-1].second;
    for (int level = 0; level <= maxLevel; level++)
    {
        for (unsigned int j = 0; j < levels.size(); j++)
        {
            // Print the identifier of the mesh
            if (level == levels[j].second)
                fprintf(lvlFile,"%d ",mesh->map_graph_elem[levels[j].first]);
        }
        fprintf(lvlFile,"\n");
    }
    fclose(lvlFile);
}

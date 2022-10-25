#ifndef CLUSTERLIST_H_
#define CLUSTERLIST_H_

#include <stdio.h>

struct cluster_t {
  int *elems;
  int nElems;
};
typedef struct cluster_t Cluster;

struct hierarchicalCluster_t {
  int branchA;
  int branchB;
  int isLeaf;
};
typedef struct hierarchicalCluster_t HierarchicalCluster;

struct integerList_t {
  int value;
  struct integerList_t *next;
};
typedef struct integerList_t IntegerList;

struct clusterList_t {
  /* Real cluster definition */
  Cluster *cl;
  int nClusters;
  int maxClusters;
  
  /* Hierarchical cluster definition */
  HierarchicalCluster *hc;
  int nHierarchicalClusters;
  IntegerList *topClusters;
};
typedef struct clusterList_t ClusterList;

ClusterList *newClusterList (int maxClusters);
void deleteClusterList (ClusterList *cl);
int newCluster (ClusterList *cl);
void addElemToCluster (ClusterList *cl, int clusterId, int elem);
int getNElemsOfCluster (ClusterList *cl, int clusterId);
int getElemOfCluster (ClusterList *cl, int clusterId, int elemId);
int getNClusters (ClusterList *cl);
int isElemInCluster (ClusterList *cl, int clusterId, int elem);
int isElemInClusterList (ClusterList *cl, int elem);
void deleteElemFromCluster (ClusterList *cl, int clusterId, int elem);
void joinClusters (ClusterList *cl, int clusterIDA, int clusterIDB);
void newHierarchicalCluster (ClusterList *cl);
void addHierarchicalCluster (ClusterList *cl, int clusterA, int clusterB);
int getNHierarchicalClusters (ClusterList *cl);
int getNElemsInTopClusters (ClusterList *cl);
int getElemInTopClusters (ClusterList *cl, int elemId);
IntegerList *getListForCluster (ClusterList *cl, int cluster);
void deleteIntegerList (IntegerList *il);
int orderInt (const void *iA, const void *iB);
void printClusterList (ClusterList *cl, FILE *fout, int firstResidue);

#endif /*CLUSTERLIST_H_*/

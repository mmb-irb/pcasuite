#include <string.h>
#include <stdlib.h>
#include "clusterlist.h"
#include "utils.h"

ClusterList *newClusterList (int maxClusters) {
  ClusterList *cl;
  
  cl = (ClusterList *)smalloc (sizeof (*cl));
  cl->maxClusters = maxClusters;
  cl->nClusters = 0;
  cl->cl = (Cluster *)smalloc (maxClusters * sizeof (*(cl->cl)));
  cl->hc = NULL;
  cl->topClusters = NULL;
  
  return cl;
}

void deleteClusterList (ClusterList *cl) {
  int i;
  
  for (i = 0; i < cl->nClusters; i++) {
    sfree (cl->cl[i].elems);
  }
  sfree (cl->cl);
  sfree (cl->hc);
  deleteIntegerList (cl->topClusters);
  sfree (cl);
}

int newCluster (ClusterList *cl) {
  cl->cl[cl->nClusters].elems = (int *)smalloc (cl->maxClusters * sizeof (int));
  cl->cl[cl->nClusters].nElems = 0;
  cl->nClusters++;
  
  return cl->nClusters - 1;
} 

void addElemToCluster (ClusterList *cl, int clusterId, int elem) {
  cl->cl[clusterId].elems[cl->cl[clusterId].nElems] = elem;
  cl->cl[clusterId].nElems++;
}

int getNElemsOfCluster (ClusterList *cl, int clusterId) {
  return cl->cl[clusterId].nElems;
}

int getElemOfCluster (ClusterList *cl, int clusterId, int elemId) {
  return cl->cl[clusterId].elems[elemId];
}

int getNClusters (ClusterList *cl) {
  if (cl == NULL)
    return 0;
  else
    return cl->nClusters;
}

int isElemInCluster (ClusterList *cl, int clusterId, int elem) {
  int found, elemId, nElems;
  
  found = FALSE;
  nElems = getNElemsOfCluster (cl, clusterId);
  elemId = 0;
  while (!found && (elemId < nElems)) {
    found = found || (getElemOfCluster (cl, clusterId, elemId) == elem);
    elemId++;
  }
  
  return found;
}

int isElemInClusterList (ClusterList *cl, int elem) {
  int found, nClusters, clusterId;
  
  found = FALSE;
  nClusters = getNClusters (cl);
  clusterId = 0;
  while (!found && (clusterId < nClusters)) {
    found = found || (isElemInCluster (cl, clusterId, elem));
    clusterId++;
  }
  
  return found;
}

void deleteElemFromCluster (ClusterList *cl, int clusterId, int elem) {
  int found, nElems, idx;
  
  /* First we search for the element in the cluster */
  nElems = getNElemsOfCluster (cl, clusterId);
  idx = 0;
  found = FALSE;
  while (!found && idx < nElems) {
    found = found || (cl->cl[clusterId].elems[idx] == elem);
    idx++;
  }

  /* Now we delete the element, shift the tail and decrement the number
   * of elements */
  if (found) {
    memmove (&(cl->cl[clusterId].elems[idx - 1]), &(cl->cl[clusterId].elems[idx]), (nElems - idx) * sizeof (int));
    cl->cl[clusterId].nElems--;
  }
}

void joinClusters (ClusterList *cl, int clusterIDA, int clusterIDB) {
  int i;
  int *clusterA, *clusterB;
  
  /* We order the clusters because we will always join the higjer cluster
   * over the lower one */
  if (clusterIDB < clusterIDA) {
    i = clusterIDB;
    clusterIDB = clusterIDA;
    clusterIDA = i;
  }
  
  /* We move the elements of one cluster onto another */
  clusterA = cl->cl[clusterIDA].elems;
  clusterB = cl->cl[clusterIDB].elems;
  for (i = 0; i < cl->cl[clusterIDB].nElems; i++) {
    clusterA[cl->cl[clusterIDA].nElems] = clusterB[i];
    cl->cl[clusterIDA].nElems++;
  }
  cl->cl[clusterIDB].nElems = 0;
  
  /* Finally, we delete the emptied cluster */
  sfree (clusterB);
  cl->cl[clusterIDB].nElems = cl->cl[cl->nClusters - 1].nElems;
  cl->cl[clusterIDB].elems = cl->cl[cl->nClusters - 1].elems;
  cl->cl[cl->nClusters - 1].nElems = 0;
  cl->nClusters--;
}

void newHierarchicalCluster (ClusterList *cl) {
  int i;
  IntegerList *il;
  
  cl->nHierarchicalClusters = 0;
  if (cl->hc != NULL) {
    sfree (cl->hc);
    deleteIntegerList (cl->topClusters);
    cl->topClusters = NULL;
  }
    
  cl->hc = (HierarchicalCluster *)smalloc (((getNClusters (cl) * 2) - 1) * sizeof (HierarchicalCluster));
  for (i = 0; i < getNClusters (cl); i++) {
    cl->hc[i].branchA = -1;
    cl->hc[i].branchB = -1;
    cl->hc[i].isLeaf = TRUE;
    il = cl->topClusters;
    cl->topClusters = (IntegerList *)smalloc (sizeof (IntegerList));
    cl->topClusters->value = i;
    cl->topClusters->next = il;
  }
  for (i = getNClusters (cl); i < ((getNClusters (cl) * 2) - 1); i++) {
    cl->hc[i].branchA = -1;
    cl->hc[i].branchB = -1;
    cl->hc[i].isLeaf = FALSE;
  }
}

void addHierarchicalCluster (ClusterList *cl, int clusterA, int clusterB) {
  int hcId;
  IntegerList *il, *aux;
  
  hcId = getNClusters (cl) + cl->nHierarchicalClusters;
  cl->hc[hcId].branchA = clusterA;
  cl->hc[hcId].branchB = clusterB;
  cl->hc[hcId].isLeaf = FALSE;
  cl->nHierarchicalClusters++;
  
  /* We add the new cluster */
  il = (IntegerList *)smalloc (sizeof (*il));
  il->value = hcId;
  il->next = cl->topClusters;
  cl->topClusters = il;
  
  /* Now we delete the old clusters (overrided by the new cluster) */
  il = cl->topClusters;
  while (il->next != NULL) {
    if (il->next->value == clusterA || il->next->value == clusterB) {
      aux = il->next;
      il->next = il->next->next;
      sfree (aux);
      continue;
    } else {
      il = il->next;
    }
  }
}

int getNElemsInTopClusters (ClusterList *cl) {
  int i;
  IntegerList *il;
  
  i = 0;
  il = cl->topClusters;
  while (il != NULL) {
    i++;
    il = il->next;
  }
  
  return i;
}

int getElemInTopClusters (ClusterList *cl, int elemId) {
  int i;
  IntegerList *il;
  
  i = 0;
  il = cl->topClusters;
  while (i != elemId) {
    i++;
    il = il->next;
  }
  
  return il->value;
}

IntegerList *getListForCluster (ClusterList *cl, int cluster) {
  IntegerList *il, *ilA, *ilB;
  int i;
  
  if (cl->hc[cluster].isLeaf) {
    il = NULL;
    for (i = 0; i < getNElemsOfCluster (cl, cluster); i++) {
      ilA = (IntegerList *)smalloc (sizeof (*il));
      ilA->value = getElemOfCluster (cl, cluster, i);
      ilA->next = il;
      il = ilA;
    } 
  } else {
    ilA = getListForCluster (cl, cl->hc[cluster].branchA);
    ilB = getListForCluster (cl, cl->hc[cluster].branchB);
    il = ilA;
    while (il->next != NULL)
      il = il->next;
    il->next = ilB;
    il = ilA;
  }
  
  return il;
}

void deleteIntegerList (IntegerList *il) {
  IntegerList *aux;

  while (il != NULL) {
    aux = il;
    il = il->next;
    sfree (aux);
  }
}

int getNHierarchicalClusters (ClusterList *cl) {
  return (getNClusters (cl) * 2) - 1;
}

int orderInt (const void *iA, const void *iB) {
  int a, b;
  
  a = *(int *)iA;
  b = *(int *)iB;

  return a - b;
}

void printClusterList (ClusterList *cl, FILE *fout, int firstResidue) {
  int i, j;
  
  /* Prior to printing, we order the lists */
  for (i = 0; i < getNClusters (cl); i++) {
    qsort (cl->cl[i].elems, getNElemsOfCluster (cl, i), sizeof (int), &orderInt);
  }
  
  /* Finally, we print the data */
  fprintf (fout, "-----ClusterList-----\n");
  fprintf (fout, "nClusters: %i\n", cl->nClusters);
  for (i = 0; i < getNClusters (cl); i++) {
    fprintf (fout, "Cluster %i (%i elements)\n  [", i, getNElemsOfCluster (cl, i));
    for (j = 0; j < getNElemsOfCluster (cl, i); j++) {
      fprintf (fout, "%i, ", getElemOfCluster (cl, i, j) + firstResidue);
    }
    fprintf (fout, "]\n");
  }
  fprintf (fout, "---------------------\n");
  
  /* VMD mask printing */
  fprintf (fout, "-----ClusterList (VMD)-----\n");
  fprintf (fout, "mol load pdb \"\"\n");
  fprintf (fout, "mol delrep 0 top\n");
  fprintf (fout, "mol color ColorId 0\n");
  fprintf (fout, "mol selection {not residue");
  for (i = 0; i < getNClusters (cl); i++) {
    for (j = 0; j < getNElemsOfCluster (cl, i); j++) {
      fprintf (fout, " %i", getElemOfCluster (cl, i, j) + firstResidue);
    }
  }
  fprintf (fout, "}\n");
  fprintf (fout, "mol addrep top\n");
  
  fprintf (fout, "mol representation Bonds\n");
  for (i = 0; i < getNClusters (cl); i++) {
    fprintf (fout, "mol color ColorId %i\n", i+1);
    fprintf (fout, "mol selection {residue");
    for (j = 0; j < getNElemsOfCluster (cl, i); j++) {
      fprintf (fout, " %i", getElemOfCluster (cl, i, j) + firstResidue);
    }
    fprintf (fout, "}\n");
    fprintf (fout, "mol addrep top\n");
  }
  fprintf (fout, "---------------------------\n");
}

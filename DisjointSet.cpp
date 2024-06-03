// DisjointSet implementation using union by size and path compression
// By Mary Elaine Califf and max glover

#include "DisjointSet.h"
#include <iostream>

DisjointSet::DisjointSet(int numObjects)
{
    numValues = numObjects;
    theArray.resize(numObjects);
    for (int i = 0; i < numObjects; i++)
    {
        theArray[i]=-1;
    }
    
}

//recursive method to find the item -- does path compression on the way out of the recursion
int DisjointSet::find(int objectIndex)
{
    int setname;
    if(theArray[objectIndex]<0){
        return objectIndex;
    }
    else{
        setname = find(theArray[objectIndex]);
        theArray[objectIndex] = setname;
        return setname;
    }
}

bool DisjointSet::doUnion(int objIndex1, int objIndex2)
{
    int val1, val2;

    val1= find(objIndex1);
    val2 = find(objIndex2);
    int numObjects = numValues*-1;

    if (theArray[val2]<theArray[val1])
    {
        theArray[val2]+=theArray[val1];
        theArray[val1]=val2;
        val1 = find(val2);
    }
    else{
        theArray[val1]+=theArray[val2];
        theArray[val2]=val1;
        val1 = find(val1);
    }

    if (numObjects== theArray[val1])
    {
     return true;   
    }
    
    return false;
}

void DisjointSet::printArrayValues(std::ostream &outputStream)
{
    for (int i = 0; i < numValues; i++)
    {
        outputStream << theArray[i] << " ";
    }
    outputStream << std::endl;
}

#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include "data2.h"

#define MAXCHILD 20 // Maxmum children the parent could have

/*! @typedef cell_t
 *  @abstract structure for each cell in the evolution history
 *  @field samid        the ID name for the celltype [char *]
 *  @field pos          the position for the celltype [int *]
 *  @field mutcount     the mutation count for the celltype [int]
*/
typedef struct _cell_t {
    char samid[60];
    int *pos;
    int mutcount;
} cell_t;

/*! @typedef cell_t
 *  @abstract structure for each cell in the evolution history
 *  @field cellnum       the cell number in the pool [int]
 *  @field cell          the cell list for the pool [cell_t *]
*/
typedef struct _pool_t {
    int cellnum;
    cell_t *cell;
} pool_t;


/*! @typedef Note
 *  @abstract structure for the evolution tree Node.
 *  @field childsize        the child count for the node [int]
 *  @field cell             the cell structure for the node
 *  @field parent           the parent pointer for the node [Node *]
 *  @field child            the children list for the node [Node **]
*/
struct _Node {
    int childsize;
    cell_t cell;
    struct _Node *parent;
    struct _Node **child;
};
typedef struct _Node Node;


/*! @funciton: creat a new matrix and allocate the necessary memory
 *   @parmeters:
 *   @    row        the row count for the new matrix [int]
 *   @    col        the col count for the new matrix [int]
 *   @return:
 *   @    the pointer to the 2D new matrix [int **]
*/
int **NewMatrix( int row, int col )
{
    int **matrix = (int **)calloc(row, sizeof(int *));
    
    if ( !matrix ) goto err;
    for ( int i=0; i < row; ++i ) {
        matrix[i] = (int *)calloc(col, sizeof(int));
        if ( !matrix[i] ) goto err;
    }
    return matrix;
    
    err:
        fprintf(stderr, "[Err::%s::%d] Faild to allocate memory!\n", __func__, __LINE__);
        return NULL;
}

/*! @funciton: free an existed matrix
 *   @parmeters:
 *   @    rawmatrix  the input matrix [int **]
 *   @    row        the row count for the new matrix [int]
 *   @    col        the col count for the new matrix [int]
 *   @return:
 *   @    void
*/
void FreeMatrix( int **matrix, int row, int col)
{
    for ( int i=0; i < row; ++i )
        free(matrix[i]);        
    free(matrix);
}


/*! @funciton: get a cell list
 *   @parmeters:
 *   @    rawmatrix  the raw input matrix [int **]
 *   @    row        the row count for the new matrix [int]
 *   @    col        the col count for the new matrix [int]
 *   @return:
 *   @    pool       a cell list [pool_t *]
*/
pool_t *GetCell( int **rawmatrix, int row, int col )
{
    pool_t *pool = (pool_t *)malloc(sizeof(pool_t)); 
    pool->cell = (cell_t *)calloc(col, sizeof(cell_t));
    
    if ( (!pool) || (!pool->cell) ) goto err;
    
    for ( int i=0; i < col; ++i ) {
        int *tempos = (int *)calloc(row, sizeof(int));
        if ( !tempos ) 
            goto err;
        else 
            (&pool->cell[i])->pos = tempos;
        
        sprintf(pool->cell[i].samid, "S%d", i+1);
        for ( int j=0; j < row; ++j ) {
            if ( rawmatrix[j][i] )
                pool->cell[i].mutcount++;
                (&pool->cell[i])->pos[j] = rawmatrix[j][i];
        }
    } pool->cellnum = col;
    return pool;
    
    err:
        fprintf(stderr, "[Err::%s::%d] Faild to allocate memory!\n", __func__, __LINE__);
        return NULL;
}

/*! @funciton: update the cell list for the pool
 *   @parmeters:
 *   @    pool        the pool for cell list [pool_t *]
 *   @    samid       the cell ID int the pool [char *]
 *   @return:
 *   @    return      0 if update sucessfully
*/
int UpdatePool( pool_t *pool, char *samid )
{
    int samindex;
    cell_t *cell = pool->cell;
    
    for ( int i=0; i < pool->cellnum; ++i ) {
        if ( !strcmp(cell[i].samid, samid) ) {
            samindex = i; break;
        }
    }
    
    pool->cellnum--; //update the cell num for the pool first

    for ( int j=samindex; j < pool->cellnum; ++j ) {
        cell[j].mutcount = cell[j+1].mutcount;
        strcpy(cell[j].samid, cell[j+1].samid);
        (&cell[j])->pos = (&cell[j+1])->pos;
    } 
    
    return 0;
}

/*! @funciton: get the cell who has the minimum snv count
 *   @parmeters:
 *   @    pool       the cell list [pool_t *]
 *   @return:
 *   @    return     the cell address who has a minimum snv count [int **]
*/
cell_t *GetMin( pool_t *pool )
{
    int minindex, mincount = 2<<25;
    cell_t *cell = pool->cell;
    
    for ( int i=0; i < pool->cellnum; ++i ) {
        if ( cell[i].mutcount < mincount ) {
            minindex = i; 
            mincount = cell[i].mutcount;
        }
    }
    
    return (&pool->cell[minindex]);
}


/*! @funciton: creat a new node and allocate the necessary memory
 *   @parmeters:
 *   @    posnum      the position num for the new node
 *   @return:
 *   @    new         the pointer to new node [Node *]
*/
Node *NewNode( void )
{
    Node *new;
    
    new = (Node *)malloc(sizeof(Node));
    new->child = (Node **)malloc(MAXCHILD * sizeof(Node*));
    if ( !new || !(new->child)) 
        goto err;
    new->cell.pos = (int *)calloc(__ROW, sizeof(int));
    if ( !new->cell.pos ) goto err;
    
    new->childsize = 0;
    new->parent = NULL;
    
    return new;
    
    err:
        fprintf(stderr, "[Err::%s::%d] Faild to allocate memory!\n", __func__, __LINE__);
        return NULL;
}

/*! @funciton: compare the snv position for the two cell
 *   @parmeters:
 *   @    T           the cell from the tree that need to be compared [cell_t *]
 *   @    cell        the cell from the pool that would be added to the tree [cell_t *]
 *   @return:
 *   @    0           the cell should have the same level with the T
 *   @    1           the cell should to be the T's child
 *   @   -1           the cell should to be the other branch
*/
int Compare( Node *point, cell_t *cell )
{
    int samnum = 0, status;
    cell_t *T = &point->cell;

    for ( int i =0; i < __ROW; ++i ) {
        if ( T->pos[i] == cell->pos[i] ) samnum++;
    }

    if ( __ROW - samnum == 0 ) // the cell is the brother of the node
        return 0;
    if ( __ROW - samnum == 2 ) {
       status = Compare( point->parent, cell); 
       if ( status == 1 ) // the cell is the brother of the node
           return 0;
       else // the cell should belong to other branch
           return -1;
    }
    if ( __ROW - samnum == 1 ) // the cell is the child of the node
        return 1;
    else // the cell should belong to other branch
        return -1;
}

/*! @funciton: search and decide which branch the new cell should be placed
 *   @parmeters:
 *   @    point       the pointer to the node [Node *]
 *   @    cell        the cell needs to be searched [cell_t *]
 *   @return:
 *   @    out         the pointer to the target node [node *]
*/
Node *Search( Node *point, cell_t *cell )
{
    Node *tem, *out, **first;
    
    if ( !point->parent && !point->childsize ) { // the root node
        return point;
    }

    if ( point->childsize == 0 ) { // only check for the leafs to speed up the search function
        int status;

        status = Compare(point, cell);
        if ( status == 0 ) // the cell should belong to 
            return (point->parent);
        else if ( status == 1 ) // the cell should belong to the child of the point
            return point;
        else // the the cell should not appear in the same branch
            return NULL;    
    }
    
    first = point->child;
    for ( int i=0; i < point->childsize; ++i ) {
        Node *tem = first[i];
        out = Search(tem, cell);
        if ( out != NULL )
            return out;
    }
    return NULL;
}

/*! @funciton: add a new node to the tree
 *   @parmeters:
 *   @    root       the pointer to the tree root node [Node *]
 *   @    cell       the cell needs to be added [cell_t *]
 *   @return:
 *   @    void
*/
void AddNode( Node *root, cell_t *cell )
{
    Node *target, *new = NewNode();
    
    new->cell.mutcount = cell->mutcount;
    strcpy(new->cell.samid, cell->samid);
    memcpy(new->cell.pos, cell->pos, __ROW * sizeof(int));
    target = Search(root, cell);
    
    new->parent = target;
    target->child[target->childsize++] = new;
        
}

/*! @funciton: display the tree infomation
 *   @parmeters:
 *   @    root        the start node of the tree [Node *]
 *   @    depth       the depth of the current node [int]
 *   @return:
 *   @    void        nothing return
*/
void PrintTree( Node *root, int depth )
{
    if ( depth == 1 ) // only display at the first node 
        printf("Root\n");

    for ( int i=0; i < depth-1; ++i ) 
        printf("|   ");

    printf("|___%s\n", root->cell.samid);

    Node **tem = root->child;
    for ( int i=0; i < root->childsize; ++i ) 
        PrintTree(tem[i], depth+1);
}


/*! @funciton: read the data from the global matrix
 *   @parmeters:
 *   @    row        the row count for the matrix [int]
 *   @    col        the col count for the matrix [int]
 *   @return:
 *   @    matrix     the raw matrix [int **]
*/
int **ReadData( int row, int col )
{
    int **matrix = NewMatrix(row, col);
    
    for ( int i=0; i < row; ++i ) {
        for ( int j=0; j < col; ++j )
            matrix[i][j] = __MATRIX[i][j];
    }
    return matrix;
}


int main( void )
{
    int **matrix, status;
    pool_t *pool;
    cell_t *curcell;
    Node *root = NewNode();

    matrix = ReadData(__ROW, __COL);
    pool = GetCell(matrix, __ROW, __COL);

    for ( int i=0; i < __COL; ++i ) {
        curcell = GetMin(pool);
        AddNode(root, curcell);
        UpdatePool(pool,curcell->samid);
    }
    printf("\n");
    PrintTree(root->child[0], 1);
    printf("\n");
}



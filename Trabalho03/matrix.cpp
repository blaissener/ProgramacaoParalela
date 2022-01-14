#include <stdlib.h>
#include <stdio.h>
#include <iostream>



int main(int argc, char *argv[]){
    
    int rows = 5, cols = 5;
    
    int** matrix = new int*[rows];
    
    for (int i = 0; i < rows; ++i){
        matrix[i] = new int[cols];
    }    

    for(int i=0; i<rows;i++){
        for(int j =0; j<cols; j++){
            matrix[i][j] = i+j;
        }
        
    }

    for(int i=0; i<rows;i++){
        for(int j =0; j<cols; j++){
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }

    for (int i = 0; i < rows; ++i)delete [] matrix[i];
    delete [] matrix;

    std::cout << "Matrix deleted" << std::endl;
    
    return 0;
}
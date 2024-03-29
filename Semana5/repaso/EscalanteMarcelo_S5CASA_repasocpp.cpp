#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include<time.h>
#include<cmath>
//#include <ctime>

// Punto 1

int x1 = 7;
double x2 = 4.3;

// Punto 8

double funcion_p8(double mivarflotante, int mivarentera){
  return mivarflotante/mivarentera;

}

// Punto 10

int min_arreglo(int * arreglo, int tamaño){
  int min = arreglo[0];
  for (int i = 0; i < tamaño; i=i+1){
    if (arreglo[i] < min){
      min = arreglo[i];
    }
  }

  return min;
}

// Punto 11
void print(int * arreglo, int tamaño){
  for (int i = 0; i<tamaño; i=i+1){
    if (arreglo[i] < 800){
      if (arreglo[i]%2 != 0){
        std::cout << arreglo[i] << " ";
      }
    }
    else {
      break;
    }
  }
}


int main(void) {

    // Punto 2
    std::cout << "La primera tiene un valor de " << x1 << " y la segunda tiene un valor de " << x2 << ".\n";
    
    // Punto 3
    double z = x2/x1;
    std::cout << "El resultado es " << z << "\n";

    // Punto 4
    const int tamaño = 300;
    int arreglo[tamaño];

    srand(time(0)); // Semilla
    for (int i = 0; i < tamaño; ++i) {
        arreglo[i] = rand() % 901;
    }

    // Punto 5
    for (int i = 0; i < tamaño; i++) {
        std::cout << arreglo[i] << " ";
    }
    std::cout << "\n";

    // Punto 6
    std::cout << arreglo[4] << "\n";

    // Punto 7
    std::cout << "La longitud del arreglo es " << sizeof(arreglo)/sizeof(arreglo[0]) << "\n";
    
    // Punto 9
    std::cout << funcion_p8(17.5, 5) << "\n";

    // Punto 10
    std::cout << min_arreglo(arreglo, tamaño) << "\n";

    // Punto 11
    print(arreglo, tamaño);

    return 0;
}


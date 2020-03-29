#include <bits/stdc++.h>
#include <ctime>
#include <chrono>
#include <random>
#define NMAX 10000002

using namespace std;
using namespace std::chrono;

ifstream fin ("teste.in");
ofstream fout("rezultate.out");

void BubbleSort(long long v[], int n)
{
    bool ok = 0;

    while (!ok) {
        ok = 1;
        for (int i = 0; i < n-1; ++i)
            if (v[i] > v[i+1]) {
                long long aux = v[i];
                v[i] = v[i+1];
                v[i+1] = aux;
                ok = 0;
            }
    }
}

void CountSort(long long v[], int n)
{
    long long val_max = v[0];
    for (int i = 1; i < n; ++i)
        if (v[i] > val_max)
            val_max = v[i];

    vector <long long> aparitii(val_max + 1, 0);
    vector <long long> output(n);

    for (int i = 0; i < n; ++i)
        aparitii[v[i]]++;

    for(int i = 1; i <= val_max; ++i)
        aparitii[i] += aparitii[i-1];

    for(int i = 0; i < n; ++i) {
        output[aparitii[v[i]]-1] = v[i];
        aparitii[v[i]]--;
    }

    for (int i = 0; i < n; ++i)
        v[i] = output[i];
}

void RadixSort(long long v[], int n)
{
    const int exp = 8; // here we can modify the exponent of the base we use (power of 2)
    int base = (1 << exp);
    vector <long long> aparitii (base);
    vector <long long> output(n);

    for (int shift = 0; shift < 64/exp; ++shift) {
        for (int i = 0; i < n; ++i)
            aparitii[(v[i] >> (exp*shift)) & (base-1)]++;

        for (int i = 1; i < base; ++i)
            aparitii[i] += aparitii[i-1];

        for (int i = n - 1; i >= 0; --i) {
            aparitii[(v[i] >> (exp*shift)) & (base-1)]--;
            output[aparitii[(v[i] >> (exp*shift)) & (base-1)]] = v[i];
        }

        for (int i = 0; i < n; ++i)
            v[i] = output[i];

        for (int i = 0; i < base; ++i)
            aparitii[i] = 0;
    }
}

void MergeFunction(long long v[], int left, int middle, int right)
{
    int lenght_one =  middle - left + 1;
    int lenght_two = right - middle;
    vector <long long> L(lenght_one);
    vector <long long> R(lenght_two);

    for(int i = 0; i < lenght_one; ++i)
        L[i] = v[left + i];

    for(int j = 0; j < lenght_two; ++j)
        R[j] = v[middle + 1 + j];

    int first_index = 0;
    int second_index = 0;
    int index = left;

    while (first_index < lenght_one && second_index < lenght_two) {
        if (L[first_index] <= R[second_index]){
            v[index] = L[first_index];
            first_index++;
        }
        else {
            v[index] = R[second_index];
            second_index++;
        }
        index++;
    }

    while (first_index < lenght_one) {
        v[index] = L[first_index];
        first_index++;
        index++;
    }

    while (second_index < lenght_two) {
        v[index] = R[second_index];
        second_index++;
        index++;
    }
}

void MergeSort(long long v[], int n, int left, int right)
{
    if (left >= right)
        return;

    int middle = left+(right-left)/2; // in case l+r is bigger than INT_SIZE
    MergeSort(v, n, left, middle);
    MergeSort(v, n, middle + 1, right);
    MergeFunction(v, left, middle, right);
}

long long MedianThree (long long v[], int left, int right){
    if (right - left + 1 < 3)
        return v[left];

    int mid = left+(right-left)/2;
    long long a = v[left];
    long long b = v[mid];
    long long c = v[right];

    if (a > b)
        swap(a,b);
    if (a > c)
        swap(a,c);
    if (b > c)
        swap(b,c);

    return b;
}

long long MedianFive (long long v[], int left, int right){

    if (right - left + 1 < 5)
        return v[left];

    long long a = v[left];
    long long b = v[left + (right - left + 1)/5];
    long long c = v[left + ((right - left + 1)/5)*2];
    long long d = v[left + ((right - left + 1)/5)*3];
    long long e = v[left + ((right - left + 1)/5)*4];

    if (a > b)
        swap(a, b);
    if (a > c)
        swap(a, c);
    if (a > d)
        swap(a, d);
    if (a > e)
        swap(a, e);
    if (b > c)
        swap(b, c);
    if (b > d)
        swap(b, d);
    if (b > e)
        swap(b, e);
    if (c > d)
        swap(c, d);
    if (c > e)
        swap (c, e);
    if (d > e)
        swap(d, e);

    return c;
}


int FixPivot(long long v[], int left, int right, long long(*MedianSelectMethod)(long long v[], int, int)){

    long long pivot = (*MedianSelectMethod)(v, left, right);
    int i;
    for (i = left; i <= right; ++i)
        if (v[i] == pivot)
            break;
    swap(v[i], v[right]);

    i = left;
    for (int j = left; j < right; ++j)
        if (v[j] <= pivot){
            swap(v[i],v[j]);
            i++;
        }
    swap(v[i], v[right]);
    return i;
}


void QuickSort (long long v[], int left, int right, long long(*MedianSelectMethod)(long long v[], int, int))
{
    if (left >= right)
        return;

    int index_pivot = FixPivot(v, left, right, (*MedianSelectMethod)); //Change the function to select the pivot (I have declared median of five aswell)

    QuickSort(v, left, index_pivot -1, (*MedianSelectMethod));
    QuickSort(v, index_pivot + 1, right, (*MedianSelectMethod));

}

void copy_vectors(long long a[], long long b[], int dim)
{
    for (int i = 0; i < dim; ++i)
        b[i] = a[i];
}


int main()
{

    int querries;
    long long *v = new long long[NMAX];
    long long *copie = new long long[NMAX];
    long long val_max, val_min;
    int dim;
    fin >> querries;
    for (int i = 1; i <= querries; ++i){
        fin >> dim >> val_min >> val_max;
        mt19937 gen(time(0));
        uniform_int_distribution<unsigned long long> dis(val_min, val_max);
        for (int j = 0; j < dim; ++j)
            v[j] = dis(gen);

        if (dim <= 10000) {
            copy_vectors(v,copie,dim);
             auto start = high_resolution_clock::now();
             BubbleSort(copie,dim);
             auto finish = high_resolution_clock::now();
             auto time = duration_cast<microseconds>(finish - start);

             if(is_sorted(copie, copie+dim))
                fout << "Sortarea folosind Bubble Sort a functionat in " << (double) time.count() << '\n';
            else
                fout << "Sortarea folosind Bubble Sort nu a functionat";

        }
        else
            fout << "Vectorul nu poate fi sortat folosind Bubble Sort deoarece ar lua prea mult timp\n";


        if (val_max <= 1000000) {
            copy_vectors(v,copie,dim);
            auto start = high_resolution_clock::now();
            CountSort(copie,dim);
            auto finish = high_resolution_clock::now();
            auto time = duration_cast<microseconds>(finish - start);

            if(is_sorted(copie,copie + dim))
                fout << "Sortarea folosind Count Sort a functionat in " << (double) time.count() << '\n';
            else
                fout << "Sortarea folosind Count Sort nu a functionat\n";
        }
        else
            fout << "Vectorul nu poate fi sortat folosind Count Sort vectorul de aparitii ar ocupa prea multa memorie\n";

        copy_vectors(v, copie,dim);
        auto start = high_resolution_clock::now();
        RadixSort(copie, dim);
        auto finish = high_resolution_clock::now();
        auto time = duration_cast<microseconds>(finish - start);

        if(is_sorted(copie,copie + dim))
            fout << "Sortarea folosind Radix Sort a functionat in " << (double) time.count() << '\n';
        else
            fout << "Sortarea folosind Radix Sort nu a functionat\n";


        copy_vectors(v,copie,dim);
        start = high_resolution_clock::now();
        MergeSort(copie, dim, 0, dim-1);
        finish = high_resolution_clock::now();
        time = duration_cast<microseconds>(finish - start);

        if(is_sorted(copie,copie + dim))
            fout << "Sortarea folosind Merge Sort a functionat in " << (double) time.count() << '\n';
        else
            fout << "Sortarea folosind Merge Sort nu a functionat\n";

        copy_vectors(v,copie,dim);
        start = high_resolution_clock::now();
        QuickSort(copie, 0, dim-1, MedianThree);
        finish = high_resolution_clock::now();
        time = duration_cast<microseconds>(finish - start);

        if(is_sorted(copie,copie + dim))
            fout << "Sortarea folosind Quick Sort(cu mediana din 3) a functionat in " << (double) time.count() << '\n';
        else
            fout << "Sortarea folosind Quick Sort(cu mediana din 3) nu a functionat\n";

        copy_vectors(v,copie,dim);
        start = high_resolution_clock::now();
        QuickSort(copie, 0, dim-1, MedianFive);
        finish = high_resolution_clock::now();
        time = duration_cast<microseconds>(finish - start);

        if(is_sorted(copie,copie + dim))
            fout << "Sortarea folosind Quick Sort(cu mediana din 5) a functionat in " << (double) time.count() << '\n';
        else
            fout << "Sortarea folosind Quick Sort(cu mediana din 5) nu a functionat\n";

        copy_vectors(v,copie,dim);
        start = high_resolution_clock::now();
        sort(copie, copie + dim);
        finish = high_resolution_clock::now();
        time = duration_cast<microseconds>(finish - start);
        fout << "Sortarea din STL a functionat in " << (double) time.count() << '\n';
        fout << '\n';
    }
    return 0;
}

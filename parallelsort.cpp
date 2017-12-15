#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <random>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <mpi.h>
#include <unistd.h>
#include <cerrno>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <cstring>
#include <iomanip>

#define NOTE_LEN 3
#define PARTS 2
#define RESULTFILE "result.txt"

int compareStr(const void* first, const void * second)
{
    const char** st1 = (const char **)first;
    const char** st2 = (const char **)second;
    return strcmp(*st1, *st2);
}

/*
void swap (  a, int* b )
{
    char* tmp;
    tmp = first;
    first = second;
    second = first;
}
*/


int partition (char** notes, int l, int h)
{
    char* x = notes[h];
    int i = l - 1;
    char* tmp;
    for (int j = l; j <= h- 1; j++)
    {
        if (compareStr(notes[j], x) < 1)
        {
            i++;
            tmp = notes[i];
            notes[i] = notes[j];
            notes[j] = tmp;
        }
    }
    tmp = notes[i+1];
    notes[i+1] = notes[h];
    notes[h] = tmp;
    return (i + 1);
}

void quickSort(char** notes, int start, int end)  //, int& wasSwap)  // if setting true > was swapped
{
    //if (start >= end)
    //    return;
    int* stack = new int[end - start + 1];
    int top = -1;
    //char* tmp;
    //int i = start, j = end;
    //char* pivot = notes[(start + end) / 2];
    
    int p;
    
    stack[++top] = start;
    stack[++top] = end;
    
    while(top >= 0)
    {
        end = stack[ top--];
        start = stack[top--];
        
        p = partition(notes, start, end);
        if(p - 1 > 1)
        {
            stack[++top] = start;
            stack[++top] = p-1;
        }
        if(p+1 < end)
        {
            stack[++top] = p + 1;
            stack[++top] = end;
        }
    }
    delete[] stack;
}

/*
        }

        while (strcmp(notes[i], pivot) < 0)
            ++i;
        while (strcmp(notes[j], pivot) > 0)
            --j;
        if (i <= j) {
            //if (i != j)
            //    wasSwap = true;*
            //std::cerr << "Swaping "<< i << " and " << j << std::endl;
            tmp = notes[i];
            notes[i] = notes[j];
            notes[j] = tmp;
            // wasSwap = true;
            // std::cerr << "Swaping "<< notes[i] << " and " << notes[j] << std::endl;
            ++i;
            --j;

           // if(i > end || j < start)
           // {
                
           // }
        }
    }
}
    
        if (j > start)
            quickSort(notes, start, j);  //, wasSwap);
        if (end > i)
            quickSort(notes, i, end);  //, wasSwap);
    //} while (i <= j);

    for(int i = start; i < end - 1; ++i)
    {
        for(int j = i + 1; j < end; ++j)
        {
            if( strcmp(notes[i],notes[j]) >= 0 )
            {
                continue;
            }

            tmp = notes[i];
            notes[i] = notes[j];
            notes[i] = tmp;
        }   

    }
}
*/

int parseArgs(int argc, char** argv, long &num, const int numProcs, std::string &filename, bool &flag)
{
    int opt;
    if (argc == 3 || argc == 5);
    else{
        std::cerr << "Invalid count of processes" << std::endl;
        return -1;
    }
    while ((opt = getopt(argc, argv, "n:f:")) != -1) {
        switch (opt) {
            case 'n':
                num = atol(optarg);
                std::cout<< "NUMBER  " << num << std::endl;
                if (num < numProcs) {
                    std::cerr << "Invalid count of notes." << std::endl;
                    return -1;
                }
                if (num % (numProcs * 2) != 0) {
                    std::cerr << "Invalid count of notes." << std::endl;
                    return -1;
                }
                break;

            case 'f':
                filename = optarg;
                flag = true;
                break;
        }
    }
    return 0;
}

int readFile(std::string filename, long &num, char** &notes, int numprocs)
{
    long count;

    std::string tmp;
    std::ifstream filestream(filename);

    if (!filestream.is_open())
    {
        std::cerr << "Open input file " << "\"" + filename + "\": " << strerror(errno) << std::endl;
        return -1;
    }
    filestream >> tmp;
    try {
        count = stol(tmp);
        if (count < numprocs) {
            std::cerr << "Error in " << "\"" + filename + "\" " <<
                "input file: " << "[1:1] "<< "Invalid count of notes." << std::endl;
            return -1;
        }
        if (count % (numprocs * 2) != 0) {
            std::cerr << "Error in " << "\"" + filename + "\" " <<
                "input file: " << "[1:1] "<< "Invalid count of notes" << std::endl;
            return -1;
        }
        num = count;
        notes = new char*[num];
        for(long i = 0; i < count; i++)
            notes[i] = new char[NOTE_LEN];

    }
    catch (std::invalid_argument&) {
        std::cerr << "Error in " << "\"" + filename + "\": " << "input file: " << "[1:]" << std::endl;
        return -1;
    }
    catch (std::bad_alloc&)
    {
        std::cerr << "Error in allocate memory" << std::endl;
        filestream.close();
        return -1;
    }

    int i = 0;
    while (filestream >> tmp) {
        if (tmp.length() != 2) {
            int position;
            if (i == 0)
                position = 2 * (i + 1);
            else
                position = 3 * (i + 1);
            std::cerr << "Error in " << "\"" + filename + "\": " << "input file: " << "[2:" << position << "]" << std::endl;
            return -1;
        }
        strcpy(notes[i], tmp.c_str());

        i ++;
    }
    if (i != count){
        std::cerr << "Error in " << "\"" + filename + "\": " << "input file: [1:]" << std::endl;
        std::cerr << "\tInvalid count of notes." << std::endl;
        std::cerr << "\tFound " << i << " notes insteed of " << count << "." << std::endl;
        return -1;
    }
    return 0;
}

int generateNotes(char** &notes, long count)
{
    srand(time(NULL));
    char* tmp = new char[NOTE_LEN];
    try {
        notes = new char*[count];
        for(long i = 0; i < count; i++)
            notes[i] = new char[NOTE_LEN];
    }
    catch (std::bad_alloc&)
    {
        std::cerr << "Error in allocate memory" << std::endl;
        return -1;
    }
    for (long i = 0; i < count; i++){
        for(int j = 0; j < NOTE_LEN; j++)
            tmp[j] = 65 + rand() % 58; // To readable ascii symbols
        tmp[NOTE_LEN - 1] = '\0';
        strcpy(notes[i], tmp);
    }
    return 0;
}

int resToFile(char** notes, long count, std::string filename)  //TODO: Change out!
{
    std::ofstream filestream(filename);
    if (!filestream.is_open())
    {
        std::cerr << "Open output file " << "\"" + filename + "\": " << strerror(errno) << std::endl;
        return -1;
    }
    filestream << count << std::endl;
    for (long i = 0; i < count; i++)
        filestream <<notes[i]<<" ";
    filestream << std::endl;
    filestream.close();
    return 0;
}

char* sendPartToArr(char** subNotes, long count, bool side)  // true -> right to array
{
    char* arr = new char[NOTE_LEN * count];
    int pos = 0;
    if (side)
    {
        for(long i = count; i < count * PARTS; i++)
            MPI_Pack(subNotes[i], NOTE_LEN, MPI_CHAR, arr, count * NOTE_LEN, &pos, MPI_COMM_WORLD);
    } 
    else
    {
        for(long i = 0; i < count; i++)
        {
            MPI_Pack(subNotes[i], NOTE_LEN, MPI_CHAR, arr, count * NOTE_LEN, &pos, MPI_COMM_WORLD);
        }
    }
    return arr;
}

void subToArr(char** subNotes, char* arr, long count)  // true -> right to array
{
    int pos = 0;
    for(long i = 0; i < count; i++)
        MPI_Pack(subNotes[i], NOTE_LEN, MPI_CHAR, arr, count * NOTE_LEN, &pos, MPI_COMM_WORLD);

}

void mergeArrays (char** subNotes, char* subArray, int subSize, bool side)  // true -> write 2 left //(sendPartToArr sends left)
{
    int pos = 0;
    char** tmp = new char*[subSize / 2];
    for (int i = 0; i < subSize / PARTS; i++)
        tmp[i] = new char[NOTE_LEN];
    if (!side) {
        for (int i = 0; i < subSize / PARTS; i++)
        {
            MPI_Unpack(subArray, subSize/2 * NOTE_LEN, &pos, tmp[i], NOTE_LEN, MPI_CHAR, MPI_COMM_WORLD);
            strcpy(subNotes[i], tmp[i]);
        }
    } else {
        for(int i = 0; i < subSize / PARTS; i++) {
            MPI_Unpack(subArray, subSize * NOTE_LEN / 2, &pos, tmp[i], NOTE_LEN, MPI_CHAR, MPI_COMM_WORLD);
            strcpy(subNotes[i + subSize/PARTS], tmp[i]);
        }
    }
}

void subShift(char** subNotes, int subSize, bool side)  //true -> shift to left
{
    if (!side)
        for (int i = subSize - 1; i >= 0; i--)
            strcpy(subNotes[i + subSize / PARTS], subNotes[i]);
    else
        for (int i = 0; i < subSize; i++)
            strcpy(subNotes[i], subNotes[subSize / PARTS + i]);
}


using namespace std;

int main (int argc, char** argv)
{
    int rank, numProcs,
            subSize,
            packPos, toRank,
            fromRank,
            maxIterations,
            swapIteration,
            currentSize;
    string filename;
    long numNotes;
    double startTime, stopTime;
    bool flag = false; // false - Exp(no file)
    
    int ErrorFromCheck = 0;

    char** notes;
    char** subNotes;
    char* array;
    char* subArray;
    char* sendArray;

    MPI_Status status;
    MPI_Request request;

    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    #define FIRST 0
    #define LAST (numProcs - 1)

    if(rank == 0) {
        ErrorFromCheck = parseArgs(argc, argv, numNotes, numProcs, filename, flag);
        if (flag && ErrorFromCheck == 0) {
            ErrorFromCheck = readFile(filename, numNotes, notes, numProcs);
        } else if (!flag && ErrorFromCheck == 0) {
            ErrorFromCheck = generateNotes(notes, numNotes);
        }
    }
    MPI_Bcast(&ErrorFromCheck, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (ErrorFromCheck != 0)
    {
        MPI_Finalize();
        return 0;
    }

    if (rank == FIRST)
    {
        //std::string inp = "input.txt";
        //resToFile(notes, numNotes, inp);
        array = new char[numNotes * NOTE_LEN];
        subToArr(notes, array, numNotes);
        subSize = numNotes / numProcs;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(&subSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    subArray = new char[subSize * NOTE_LEN];  //TODO TRY

    subNotes = new char*[int(subSize * 1.5)];
    for(int i = 0; i < int(subSize * 1.5); i++)
        subNotes[i] = new char[NOTE_LEN];

    currentSize = subSize;

    MPI_Scatter(array, subSize * NOTE_LEN, MPI_PACKED, subArray, subSize * NOTE_LEN, MPI_PACKED, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    packPos = 0;
    for(int i = 0; i < subSize; i++)
        MPI_Unpack(subArray, subSize * NOTE_LEN, &packPos, subNotes[i], NOTE_LEN, MPI_CHAR, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == FIRST)
        startTime = MPI_Wtime();

    swapIteration = 1;  //send to next
    maxIterations = (numProcs % 2 == 0) ? (numProcs + 1) : (numProcs + 2);
    do {
        bool side;

        std::qsort(subNotes, currentSize, sizeof(char*), compareStr);
        //quickSort(subNotes, 0, currentSize - 1);

        if (swapIteration % 2 == 0)
        {
            side = true;  //  send to left!
            if (FIRST == LAST)
            {

            }
            else if (rank == FIRST)
            {   
                fromRank = FIRST + 1;
                sendArray = new char[NOTE_LEN * subSize / 2];
                currentSize = subSize;
                MPI_Irecv(sendArray, (subSize / 2) * NOTE_LEN, MPI_CHAR, fromRank, rank + fromRank, MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);
                mergeArrays(subNotes, sendArray, subSize, side);
            }
            else if (rank == LAST)
            {
                toRank = LAST - 1;
                sendArray = sendPartToArr(subNotes, subSize / PARTS, !side);
                MPI_Isend(sendArray, (subSize / PARTS) * NOTE_LEN, MPI_CHAR, toRank, rank + toRank,
                          MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);
                subShift(subNotes, subSize, side);
                currentSize = subSize;
            }
            else if (rank != LAST && rank != FIRST)
            {
                toRank = rank - 1;
                fromRank = rank + 1;
                sendArray = sendPartToArr(subNotes, subSize / PARTS, !side);
                subShift(subNotes, subSize / 2, side);
                MPI_Isend(sendArray, (subSize / PARTS) * NOTE_LEN, MPI_CHAR, toRank, rank + toRank,
                          MPI_COMM_WORLD, &request);
                MPI_Irecv(sendArray, (subSize / PARTS) * NOTE_LEN, MPI_CHAR, fromRank, rank + fromRank,
                          MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);
                mergeArrays(subNotes, sendArray, subSize, side);
            }
        }
        else
        {
            side = false;  // send to right

            if (FIRST == LAST)
            {

            }
            else if (rank == FIRST)
            {
                toRank = FIRST + 1;
                sendArray = sendPartToArr(subNotes, subSize / PARTS, !side);
                currentSize = subSize / PARTS;

                MPI_Isend(sendArray, currentSize * NOTE_LEN, MPI_CHAR, toRank, rank + toRank,
                      MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);
            }
            else if (rank == LAST)
            {
                fromRank = LAST - 1;
                sendArray = new char[subSize / 2 * NOTE_LEN];
                subShift(subNotes, subSize, side);
                currentSize = subSize * 1.5;
                MPI_Irecv(sendArray, (subSize / 2) * NOTE_LEN, MPI_CHAR, fromRank, rank + fromRank, MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);
                mergeArrays(subNotes, sendArray, subSize, side);
            }
            else if (rank != LAST && rank != FIRST)
            {
                toRank = rank + 1;
                fromRank = rank - 1;
                sendArray = sendPartToArr(subNotes, subSize / PARTS, !side);
                subShift(subNotes, subSize/PARTS, side);
                MPI_Isend(sendArray, subSize / 2 * NOTE_LEN, MPI_CHAR, toRank, rank + toRank,
                          MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);
                MPI_Irecv(sendArray, (subSize / PARTS) * NOTE_LEN, MPI_CHAR, fromRank, rank + fromRank,
                          MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);
                mergeArrays(subNotes, sendArray, subSize, side);
            }

        }
        
        //quickSort(subNotes, 0, currentSize -1);  //, wasSwaped);

        ++swapIteration;
    } while (swapIteration < maxIterations);
    
    
    if(rank == FIRST)
        stopTime = MPI_Wtime();

    subToArr(subNotes, subArray, subSize);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(subArray, subSize * NOTE_LEN, MPI_CHAR, array, subSize * NOTE_LEN, MPI_CHAR, FIRST, MPI_COMM_WORLD);

    if(rank == FIRST)
    {
        int pos = 0;
        for (int i = 0; i < numNotes; i++)
            MPI_Unpack(array, numNotes * NOTE_LEN, &pos, notes[i], NOTE_LEN, MPI_CHAR, MPI_COMM_WORLD);

        if(flag)
            resToFile(notes, numNotes, RESULTFILE);
        else
	{
            cout << "Data: " << numNotes << " notes. Procs: " << numProcs << "." << endl;
            cout << "Elapsed time: " << std::fixed << std::setprecision(2) << stopTime - startTime << endl;
            //resToFile(notes, numNotes, RESULTFILE);
        }
        delete[] array;
        for(int i = 0; i < numNotes; i++)
            delete[] notes[i];
        delete[] notes;
    }

    if(FIRST != LAST)
    {
        delete[] sendArray;
    }
    
    for(int i = 0; i < subSize + subSize / PARTS; i++)
	delete[] subNotes[i];
    delete[] subNotes;

    MPI_Finalize();

    return 0;
}

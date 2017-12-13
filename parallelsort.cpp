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

#define NOTE_LEN 3
#define PARTS 2
#define RESULTFILE "result.txt"

double startT, stopT;

void quickSort(char** notes, int start, int end)  //, int& wasSwap)  // if setting true > was swapped
{
    if (start >= end)
        return;
    char* tmp;
    int i = start, j = end;
    char* pivot = notes[(start + end) / 2];
    do {
        while (strcmp(notes[i], pivot) < 0)
            ++i;
        while (strcmp(notes[j], pivot) > 0)
            --j;
        if (i <= j) {
            /*if (i != j)
                wasSwap = true;*/
            //std::cerr << "Swaping "<< i << " and " << j << std::endl;
            tmp = notes[i];
            notes[i] = notes[j];
            notes[j] = tmp;
            // wasSwap = true;
            // std::cerr << "Swaping "<< notes[i] << " and " << notes[j] << std::endl;
            ++i;
            --j;
        }
        if (j > start)
            quickSort(notes, start, j);  //, wasSwap);
        if (end > i)
            quickSort(notes, i, end);  //, wasSwap);
    } while (i <= j);
}

int parseArgs(int argc, char** argv, int &num, const int numProcs, std::string &filename, bool &flag)
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
                num = atoi(optarg);
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

int readFile(std::string filename, int &num, char** &notes, int numprocs)
{
    int count;
//    bool flag = true;  // true if count % numprocs != 0
    std::string tmp;
    std::ifstream filestream(filename);

    if (!filestream.is_open())
    {
        std::cerr << "Open input file " << "\"" + filename + "\": " << strerror(errno) << std::endl;
        return -1;
    }
    filestream >> tmp;
    try {
        count = stoi(tmp);
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
        for(int i = 0; i < count; i++)
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
        std::cerr << "Error in " << "\"" + filename + "\": " << "input file:" << std::endl;
        std::cerr << "\tInvalid count of notes." << std::endl;
        std::cerr << "\tFound " << i << " notes insteed of " << count << "." << std::endl;
        return -1;
    }
//    for (int i = 0; i < num; i++)
//        std::cout<<" "<<"SubNote "<<i<<":"<<notes[i]<<std::endl;
    return 0;
}

int generateNotes(char** &notes, int count)
{
    srand(time(NULL));
    char* tmp = new char[NOTE_LEN];
    try {
        notes = new char*[count];
        for(int i = 0; i < count; i++)
            notes[i] = new char[NOTE_LEN];
    }
    catch (std::bad_alloc&)
    {
        std::cerr << "Error in allocate memory" << std::endl;
        return -1;
    }
    for (int i = 0; i < count; i++){
        for(int j = 0; j < NOTE_LEN; j++)
            tmp[j] = 65 + rand() % 58; //  rand() % 256;  //65 + rand() % 58;
        tmp[NOTE_LEN - 1] = '\0';
        strcpy(notes[i], tmp);
        //std::cout<<"NOTE "<<i<<":"<<notes[i]<<std::endl;
    }
    return 0;
}

int resToFile(char** notes, int count, std::string filename)  //TODO: Change out!
{
    std::ofstream filestream(filename);
    if (!filestream.is_open())
    {
        std::cerr << "Open output file " << "\"" + filename + "\": " << strerror(errno) << std::endl;
        return -1;
    }
    filestream << count << std::endl;
    for (int i = 0; i < count; i++)
        filestream <<notes[i]<<" ";
    filestream.close();
    return 0;
}

char* sendPartToArr(char** subNotes, int count, bool side)  // true -> right to array
{
//        std::cout << "LEN IS" << count << ":" << NOTE_LEN * count << std::endl;

    char* arr = new char[NOTE_LEN * count];
    int pos = 0;
    if (side) {
        for(int i = count; i < count * PARTS; i++)
            MPI_Pack(subNotes[i], NOTE_LEN, MPI_CHAR, arr, count * NOTE_LEN, &pos, MPI_COMM_WORLD);
    } else {
        for(int i = 0; i < count; i++) {
            //std::cout<<"PACKING: "<<"NOTE "<<i<<":"<<subNotes[i]<<std::endl;
            MPI_Pack(subNotes[i], NOTE_LEN, MPI_CHAR, arr, count * NOTE_LEN, &pos, MPI_COMM_WORLD);
        }
            //MPI_Pack(subNotes[i], NOTE_LEN, MPI_CHAR, arr, count * NOTE_LEN, &pos, MPI_COMM_WORLD);
    }
    return arr;
}

void subToArr(char** subNotes, char* arr, int count)  // true -> right to array
{
    int pos = 0;
    for(int i = 0; i < count; i++)
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
            //std::cout<<" "<<"ITARATE! "<<i<<":"<<subNotes[i]<<std::endl;
            MPI_Unpack(subArray, subSize/2 * NOTE_LEN, &pos, tmp[i], NOTE_LEN, MPI_CHAR, MPI_COMM_WORLD);
            strcpy(subNotes[i], tmp[i]);
        }
    } else {
        for(int i = 0; i < subSize / PARTS; i++) {
            MPI_Unpack(subArray, subSize * NOTE_LEN / 2, &pos, tmp[i], NOTE_LEN, MPI_CHAR, MPI_COMM_WORLD);
            strcpy(subNotes[i + subSize/PARTS], tmp[i]);
        }
    }
    //delete []tmp;
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
            numNotes, subSize,
            packPos, toRank,
            fromRank,
            maxIterations,
            swapIteration,
            currentSize;
    string filename;
    bool flag = false; // false - Exp(no file)
    //int wasSwaped;
    //int wasSwappedCheck;
    //bool needSort;
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
            //std::cout<<"Note "<<i<<":"<<notes[i]<<std::endl;
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
    std::string inp = "input.txt";
    resToFile(notes, numNotes, inp);

    if(rank == 0) {
        array = new char[numNotes * NOTE_LEN];
        currentSize = subSize;

//        for (int i = 0; i < numNotes; i++)
//            std::cout<<" "<<"SubNote "<<i<<":"<<notes[i]<<std::endl;

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

    //needSort = 0;
    MPI_Scatter(array, subSize * NOTE_LEN, MPI_PACKED, subArray, subSize * NOTE_LEN, MPI_PACKED, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    packPos = 0;
    for(int i = 0; i < subSize; i++)
        MPI_Unpack(subArray, subSize * NOTE_LEN, &packPos, subNotes[i], NOTE_LEN, MPI_CHAR, MPI_COMM_WORLD);

    quickSort(subNotes, 0, subSize-1); //, wasSwaped);

    swapIteration = 1;  //send to next
    maxIterations = (numProcs % 2 == 0) ? (numProcs + 3) : (numProcs + 2);
    do {
        bool side;

        /*if (rank == LAST)
        {
            //std::cout<<"ITER:"<<swapIteration<<std::endl;
        }*/

        if (swapIteration % 2 == 0)
        {
            side = true;  //  send to left!
            if (FIRST == LAST)
            {

            }
            else if (rank == FIRST)
            {
                //int tag;
                if(numProcs > 1)
                    fromRank = FIRST + 1;
                else
                    fromRank = FIRST;
                sendArray = new char[NOTE_LEN * subSize / 2];
                /*for (int i = 0; i < currentSize; i++)
                    std::cout<<"BEFORE RECV "<<swapIteration<<" :"<<rank<<" "<<"NOTE "<<i<<":"<<subNotes[i]<<std::endl;*/
                currentSize = subSize;
                MPI_Irecv(sendArray, (subSize / 2) * NOTE_LEN, MPI_CHAR, fromRank, rank + fromRank,
                          MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);
                mergeArrays(subNotes, sendArray, subSize, side);
            }
            else if (rank == LAST)
            {
                if(numProcs > 1)
                    toRank = LAST - 1;
                else
                    toRank = LAST;
                //int tag = rank + toRank;
                sendArray = sendPartToArr(subNotes, subSize / PARTS, !side);
                MPI_Isend(sendArray, (subSize / PARTS) * NOTE_LEN, MPI_CHAR, toRank, rank + toRank,
                          MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);
                /*for (int i = 0; i < currentSize; i++)
                    std::cout<<"BEFORE shift:"<<rank<<" "<<"NOTE "<<i<<":"<<subNotes[i]<<std::endl;*/
                subShift(subNotes, subSize, side);
                currentSize = subSize;
                /*for (int i = 0; i < currentSize; i++)
                    std::cout<<"AFTER SHIRT:"<<rank<<" "<<"NOTE "<<i<<":"<<subNotes[i]<<std::endl;*/
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
        else  //std::cout<<"rank:"<<rank<<" "<<"SubNote "<<sendArray<<std::endl;
        {
            side = false;  // send to right

            if (FIRST == LAST)
            {

            }
            else if (rank == FIRST)
            {
                if(numProcs > 1)
                    toRank = FIRST + 1;
                else
                    toRank = FIRST;
//                toRank = FIRST + 1;
                /*for (int i = 0; i < currentSize; i++)
                    std::cout<<"BEFORE SEND "<<swapIteration<<" :"<<rank<<" "<<"NOTE "<<i<<":"<<subNotes[i]<<std::endl;*/
                sendArray = sendPartToArr(subNotes, subSize / PARTS, !side);
                currentSize = subSize / PARTS;

                MPI_Isend(sendArray, currentSize * NOTE_LEN, MPI_CHAR, toRank, rank + toRank,
                      MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);
                //std::cout<<"rank:"<<rank<<" "<<"I AM  "<<rank<<std::endl;
                /*for (int i = 0; i < currentSize; i++)
                    std::cout<<"AFTER SEND "<<swapIteration<<" :"<<rank<<" "<<"NOTE "<<i<<":"<<subNotes[i]<<std::endl;*/
            }
            else if (rank == LAST)
            {
                if(numProcs > 1)
                    fromRank = LAST - 1;
                else
                    fromRank = LAST;
                //fromRank = LAST - 1;
                sendArray = new char[subSize / 2 * NOTE_LEN];
                /*for (int i = 0; i < currentSize; i++)
                    std::cout<<"BEFORE SHIFT:"<<rank<<" "<<"NOTE "<<i<<":"<<subNotes[i]<<std::endl;*/
                subShift(subNotes, subSize, side);
                currentSize = subSize * 1.5;
                /*for (int i = 0; i < currentSize; i++)
                    std::cout<<"AFTER SHIFT:"<<rank<<" "<<"NOTE "<<i<<":"<<subNotes[i]<<std::endl;*/
                MPI_Irecv(sendArray, (subSize / 2) * NOTE_LEN, MPI_CHAR, fromRank, rank + fromRank,
                          MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);

                /*for (int i = 0; i < currentSize; i++)
                    std::cout<<"OLD:"<<rank<<" "<<"SubNote "<<i<<":"<<subNotes[i]<<std::endl;*/
                mergeArrays(subNotes, sendArray, subSize, side);
                /*for (int i = 0; i < currentSize; i++)
                    std::cout<<"AFTER MERGE:"<<rank<<" "<<"NOTE "<<i<<":"<<subNotes[i]<<std::endl;
                for (int i = 0; i < currentSize; i++)
                    std::cout << "NEW:" << i << " " << "SubNote " << i << ":" << subNotes[i] << std::endl;*/
                //std::cout<<"rank:"<<rank<<" "<<"I AM  "<<rank<<std::endl;
            }
            else if (rank != LAST && rank != FIRST)
            {
                toRank = rank + 1;
                fromRank = rank - 1;
//                for (int i = 0; i < subSize; i++)
//                    std::cout << "NEW:" << i << " " << "SubNote " << i << ":" << subNotes[i] << std::endl;
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

        quickSort(subNotes, 0, currentSize -1);  //, wasSwaped);

/*        MPI_Reduce(&wasSwaped, &wasSwappedCheck, 1, MPI_INT, MPI_SUM, FIRST, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == FIRST) {
            //std::cout<<"Reduce OP: "<<wasSwappedCheck<<std::endl;
            if (wasSwappedCheck != 0)
                needSort = true;
            else
                needSort = false;
            wasSwappedCheck = 0;
        }
        MPI_Bcast(&needSort, 1, MPI_BYTE, 0, MPI_COMM_WORLD);

        for (int i = 0; i < subSize; i++)
            std::cout<<"rank:"<<rank<<" "<<"NOTE "<<i<<":"<<subNotes[i]<<std::endl;

        MPI_Finalize();
        return 0;*/
        /*if (rank == LAST)
        {
            std::cout<<"ITER:"<<swapIteration<<std::endl;
        }*/
            //std::cout<<"ITER:"<<swapIteration<<std::endl;
        ++swapIteration;
    } while (swapIteration < maxIterations);

//    for (int i = 0; i < currentSize; i++)
//        std::cout<<"NEW:"<<rank<<" "<<"SubNote "<<i<<":"<<subNotes[i]<<std::endl;

    subToArr(subNotes, subArray, subSize);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(subArray, subSize * NOTE_LEN, MPI_CHAR, array, subSize * NOTE_LEN, MPI_CHAR, FIRST, MPI_COMM_WORLD);

    //std::cout<<"HUYAK:"<<rank<<" "<<std::endl;

    int pos = 0;
    for (int i = 0; i < numNotes; i++)
        MPI_Unpack(array, numNotes * NOTE_LEN, &pos, notes[i], NOTE_LEN, MPI_CHAR, MPI_COMM_WORLD);


    if (rank == FIRST)
        resToFile(notes, numNotes, RESULTFILE);
    MPI_Finalize();
    /*
        startT = MPI_Wtime(); //Start The Time.
        stopT = MPI_Wtime();
        double parallelTime = stopT- startT;
        printf("\n\n\nTime: %f", parallelTime);

        double poslTime = stopT - startT;
        printf("\n Time: %f", stopT - startT);
        printf("\n SpeedUp: %f \n\n\n", poslTime/parallelTime);

     */
    return 0;
}
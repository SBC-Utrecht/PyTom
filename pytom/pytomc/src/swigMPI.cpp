#include "swigMPI.hpp"

#include <stdlib.h>
#include <mpi.h>
#include <stdexcept>
#include <string>
#include <memory>
#include <sstream>

void init(){

	int error = MPI_Init(NULL,NULL);
	if (error != MPI_SUCCESS) {
		throw std::runtime_error("Error initialising MPI. ABORTING!");
	}
}

bool isInitialised(){
	int flag;
	int error = MPI_Initialized(&flag);

	if (error != MPI_SUCCESS) {
		throw std::runtime_error("Could not determine whether MPI is running!");
	}

	return flag;

}

int size(){

	int mpiSize;
	int error = MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

	if (error != MPI_SUCCESS) {
			throw std::runtime_error("Could not determine COMM size!");
		}

	return mpiSize;
}

int rank(){

	int id=-1,error;

	error = MPI_Comm_rank(MPI_COMM_WORLD,&id);
	if (error != MPI_SUCCESS) {
		throw std::runtime_error("Could not determine machine rank!");
	}
	return id;
}

void send(std::string message,int recipient){

	int messageSize = message.length();
	//send data how big the string is
	int error = MPI_Send((void*)&messageSize,1,MPI_INT,recipient,0,MPI_COMM_WORLD);

	if (error != MPI_SUCCESS) {
		throw std::runtime_error("Error sending metadata message!");
	}
	//send real string
	error = MPI_Send((void *) message.c_str(),  messageSize , MPI_CHAR, recipient,0, MPI_COMM_WORLD);

	if (error != MPI_SUCCESS) {
		throw std::runtime_error("Error sending message!");
	}
}

std::string receive(){
	int messageSize=0;
	MPI_Status status,status2;
	//determine how big the string is
	int error = MPI_Recv((void *)&messageSize, 1, MPI_INT,MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

	if (error != MPI_SUCCESS) {
		throw std::runtime_error("Error receiving metadata message!");
	}

	std::auto_ptr<char> buffer((char*)malloc(sizeof(char)*messageSize));
	//receive the string
	error = MPI_Recv((void*)buffer.get(), messageSize, MPI_CHAR,status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status2);

	if(!status.MPI_SOURCE == status2.MPI_SOURCE){
		throw std::runtime_error("Concurency error in messages. Received main message from different sender than metadata message!");
	}
	if (error != MPI_SUCCESS) {
		throw std::runtime_error("Error receiving message!");
	}

	return std::string(buffer.get(),messageSize);
}

void finalise(){

	int error = MPI_Finalize();
	if (error != MPI_SUCCESS){
		throw std::runtime_error("Error finalizing MPI environment!");
	}


}

/*
 * swigMPI.hpp
 *
 *  Created on: Feb 10, 2010
 *      Author: Thomas Hrabe
 */
#ifndef SWIGMPI_HPP_
#define SWIGMPI_HPP_

#include <string>
void init();

bool isInitialised();

int rank();
int size();
void send(std::string message,int recipient);

std::string receive();
//void Comm::Recv(void* buf, int count, const Datatype& datatype,int source, int tag, Status& status) const

void finalise();

#endif /* SWIGMPI_HPP_ */

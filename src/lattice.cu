
#include "../include/lattice.h"
#include <iostream>

using namespace LatticeGraph;

// Constructor and destructor

Lattice::Lattice(){
}

Lattice::~Lattice(){
	this->getVortices().clear();
	this->getEdges().clear();
}

// Getters

std::vector< std::shared_ptr<Node> >& Lattice::getVortices(){
	return this->vortices;
}

std::shared_ptr<Node> Lattice::getVortexIdx(unsigned int idx){
	return getVortices().at(idx);
}

/***
 * Gets the location of the Node with UID uid.
 */
unsigned int Lattice::getVortexIdxUid(unsigned int uid){
	for (size_t ii=0; ii< getVortices().size(); ++ii){
		if(this->Lattice::getVortexIdx(ii)->getUid()== uid){
			return ii;
		}
	}
	return -1;
}

/***
 * Gets the the Node with UID uid. Assumes Node exists.
 */
std::shared_ptr<Node> Lattice::getVortexUid(unsigned int uid){
	for (std::shared_ptr<Node> n : this->Lattice::getVortices()){
		if(n->getUid()== uid){
			return n;
		}
	}
	return std::shared_ptr<Node>();
}

double Lattice::getVortexDistance(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2){
	return sqrt(pow(n1->getData().getCoords().x - n2->getData().getCoords().x,2)
	            +  pow(n1->getData().getCoords().y - n2->getData().getCoords().y,2));
}

double Lattice::getVortexDistanceD(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2){
	return sqrt(pow(n1->getData().getCoordsD().x - n2->getData().getCoordsD().x,2)
	            +  pow(n1->getData().getCoordsD().y - n2->getData().getCoordsD().y,2));
}

std::shared_ptr<Edge> Lattice::getEdgeIdx(unsigned int idx){
	return getEdges().at(idx);
}

/***
 * Gets the location of the Edge with UID uid.
 */
unsigned int Lattice::getEdgeIdxUid(unsigned int uid){
	for (size_t ii=0; ii< getEdges().size(); ++ii){
		if(this->Lattice::getEdgeIdx(ii)->getUid()== uid){
			return ii;
		}
	}
	return -1;
}

/***
 * Gets the the Edge with UID uid. Assumes Node exists.
 */
std::shared_ptr<Edge> Lattice::getEdgeUid(unsigned int uid){
	for (std::shared_ptr<Edge> e : this->Lattice::getEdges()){
		if(e->getUid()== uid){
			return e;
		}
	}
	return NULL;
}

std::vector< std::shared_ptr<Edge> >& Lattice::getEdges(){
	return this->edges;
}

// Setters

void Lattice::setVortex(unsigned int idx, std::shared_ptr<Node> n){
	this->Lattice::getVortices().at(idx)=(n);
}

void Lattice::setEdge(unsigned int idx, std::shared_ptr<Edge> e){
	this->Lattice::getEdges().at(idx)=(e);
}

// Creation


void Lattice::createEdges(unsigned int radius){
	std::shared_ptr<Edge> e;
	double dist = 0.0;
	for(size_t ii=0; ii< this->Lattice::getVortices().size(); ++ii){
		//std::cout << "Got here ii " << ii << std::endl;
		for(size_t jj=ii+1; jj < this->Lattice::getVortices().size(); ++jj){
			dist = Lattice::getVortexDistance(this->getVortexIdx(ii),this->getVortexIdx(jj));
			if(dist < radius ) {
				//std::cout << "Got here jj " << jj << std::endl;
				e.reset(new Edge ( this->getVortexIdx(ii), this->getVortexIdx(jj) ));
				e->setWeight(dist);
				this->Lattice::addEdge(e,this->getVortexIdx(ii),this->getVortexIdx(jj));
			}
		}
	}
}
void Lattice::createEdges(double radius){
	std::shared_ptr<Edge> e;
	double dist = 0.0;
	for(size_t ii=0; ii< this->Lattice::getVortices().size(); ++ii){
		//std::cout << "Got here ii " << ii << std::endl;
		for(size_t jj=ii+1; jj < this->Lattice::getVortices().size(); ++jj){
			dist = Lattice::getVortexDistance(this->getVortexIdx(ii),this->getVortexIdx(jj));
			if( dist < radius ) {
				//std::cout << "Got here jj " << jj << std::endl;
				e.reset(new Edge ( this->getVortexIdx(ii), this->getVortexIdx(jj) ));
				e->setWeight(dist);
				this->Lattice::addEdge(e,this->getVortexIdx(ii),this->getVortexIdx(jj));
			}
		}
	}
}

void Lattice::addVortex(std::shared_ptr<Node> n){
	this->Lattice::getVortices().push_back((n));
}

void Lattice::addEdge(std::shared_ptr<Edge> e){
	this->addEdge(e, e->getVortex(0).lock(), e->getVortex(1).lock());
}

void Lattice::addEdge(std::shared_ptr<Edge> e, std::shared_ptr<Node> n1, std::shared_ptr<Node> n2){
	this->Lattice::getEdges().push_back(e);
	std::weak_ptr<Edge> e1 = e;
	std::weak_ptr<Edge> e2 = e;
	n1->addEdge(e1);
	n2->addEdge(e2);
}

// Deletion

void Lattice::removeVortex(std::shared_ptr<Node> n){
	for(std::weak_ptr<Edge> e : n->getEdges()){
		if(e.lock()){
			std::cout << "UID: Removing Vortex{" << n->getUid() <<"}." << std::endl;
			this->removeEdge(e.lock());
			this->Lattice::getVortices().erase(this->Lattice::getVortices().begin() + this->getVortexIdxUid(n->getUid()));
		}
		else{
			std::cout << "Cannot remove UID:Edge{"<< e.lock()->getUid() << "}, does not exist." << std::endl;
		}
	}
}

void Lattice::removeVortexUid(unsigned int uid){
	auto vtx = this->getVortexUid(uid);
	if(vtx){
		this->Lattice::removeVortex(vtx);
	}
	else{
		std::cout << "Cannot remove UID:Vortex{"<< uid << "}, does not exist." << std::endl;
	}
}

void Lattice::removeVortexIdx(unsigned int idx){
	auto vtx = this->getVortexIdx(idx);
	if(vtx){
		this->Lattice::removeVortex(vtx);
	}
	else{
		std::cout << "Cannot remove IDX:Vortex["<< idx << "], does not exist." << std::endl;
	}
}

void Lattice::removeEdge(std::shared_ptr<Edge> e){
	std::cout << "Removing Edge{" << e->getUid() <<"} connecting Node{" << e->getVortex(0).lock()->getUid() << "} and Node{" << e->getVortex(1).lock()->getUid() << "}." << std::endl;
	e->getVortex(0).lock()->removeEdgeUid(e->getUid());
	e->getVortex(1).lock()->removeEdgeUid(e->getUid());
	this->Lattice::getEdges().erase(this->Lattice::getEdges().begin() + this->Lattice::getEdgeIdxUid(e->getUid()));
}

void Lattice::removeEdgeIdx(unsigned int idx){
	std::weak_ptr<Edge> e = this->getEdgeIdx(idx);
	if (auto el = e.lock()) {
		this->Lattice::removeEdge(el);
	}
	else{
		std::cout << "Cannot remove IDX:Edge[" << idx << "], does not exist." << std::endl;
	}
}

void Lattice::removeEdgeUid(unsigned int uid) {
	std::weak_ptr<Edge> e = this->getEdgeUid(uid);
	if (auto el = e.lock()) {
		this->Lattice::removeEdge(el);
	}
	else{
		std::cout << "Cannot remove UID:Edge{" << uid << "}, does not exist." << std::endl;
	}
}

void Lattice::removeEdge(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2){
	std::weak_ptr<Edge> e = this->Lattice::isConnected(n1,n2);
	if(e.lock()){
		this->Lattice::removeEdge(e.lock());
	}
	else{
		std::cout << "Node{" << n1->getUid() << "} and Node{" << n2->getUid() << "} were unconnected." << std::endl;
	}

}

void Lattice::removeEdges(std::shared_ptr<Node> n1){
	//n1->removeEdges();
}


void Lattice::createVortex(double posx, double posy, int winding){

}

void Lattice::destroyVortex(unsigned int uid){
	this->Lattice::getVortexUid(uid);
}

// Generating matrices

/**
 * Create adjacency matrix
 */
void Lattice::genAdjMat(unsigned int *mat){
	int idx1, idx2, idx;
	idx1 = 0; idx2 = 0; idx=0;
	for(std::shared_ptr<Node> n : this->Lattice::getVortices()){
		idx1=this->getVortexIdxUid(n->getUid());
		for(std::weak_ptr<Edge> e : n->getEdges()){
			idx2 = this->getVortexIdxUid(n->getConnectedNode(e.lock())->getUid());
			//std::cout << "this=" << n->getUid() << "   connected=" << n->getConnectedNode(e.lock())->getUid() << std::endl;
			idx = idx1*this->Lattice::getVortices().size() + idx2;
			//std::cout << "idx1=" << idx1 << "   idx2=" << idx2 << " idx=" << idx << "\n" << std::endl;
			mat[idx] = 1;
		}
	}
}

void Lattice::genAdjMat(double *mat){
	int idx1, idx2, idx;
	idx1 = 0; idx2 = 0; idx=0;
	for(std::shared_ptr<Node> n : this->Lattice::getVortices()){
		idx1=this->getVortexIdxUid(n->getUid());
		for(std::weak_ptr<Edge> e : n->getEdges()){
			idx2 = this->getVortexIdxUid(n->getConnectedNode(e.lock())->getUid());
			//std::cout << "this=" << n->getUid() << "   connected=" << n->getConnectedNode(e.lock())->getUid() << std::endl;
			idx = idx1*this->Lattice::getVortices().size() + idx2;
			//std::cout << "idx1=" << idx1 << "   idx2=" << idx2 << " idx=" << idx << "\n" << std::endl;
			mat[idx] = this->Lattice::getVortexDistance(n, this->getVortexIdx(idx2));
		}
	}
}

/**
 * Outputs adjacency matrix in format for copy/paste into Mathematica.
 */
void Lattice::adjMatMtca(unsigned int *mat){
	unsigned int size = this->Lattice::getVortices().size();
	std::cout << "{";
	for(size_t ii = 0; ii < size; ++ii){
		std::cout << "{";
		for(size_t jj = 0; jj < size; ++jj){
			std::cout << mat[ii*size + jj];
			if(jj<size-1)
				std::cout <<",";
			else
				std::cout << "}";
		}
		if(ii<size-1)
			std::cout <<",";
		std::cout << std::endl;
	}
	std::cout << "}" << std::endl;
}
void Lattice::adjMatMtca(double *mat){
	unsigned int size = this->Lattice::getVortices().size();
	std::cout << "{";
	for(size_t ii = 0; ii < size; ++ii){
		std::cout << "{";
		for(size_t jj = 0; jj < size; ++jj){
			std::cout << mat[ii*size + jj];
			if(jj<size-1)
				std::cout <<",";
			else
				std::cout << "}";
		}
		if(ii<size-1)
			std::cout <<",";
		std::cout << std::endl;
	}
	std::cout << "}" << std::endl;
}

// Check connection

std::weak_ptr<Edge> Lattice::isConnected(std::shared_ptr<Node> n1, std::shared_ptr<Node> n2){

	if(n1->getUid() != n2->getUid()){
		for(std::weak_ptr<Edge> e1 : n1->getEdges()){
			if(e1.lock()->isMember(n2)){
				return e1;
			}
		}
	}
	return std::weak_ptr<Edge> ();
}

// Swapping indices and ids

void Lattice::swapIdxUid(unsigned int uid1, unsigned int uid2) {
	Lattice::swapIdx(this->getVortexIdxUid(uid1),this->getVortexIdxUid(uid2));
}
void Lattice::swapIdx(unsigned int idx1, unsigned int idx2) {
	std::swap(this->getVortices().at(idx1),this->getVortices().at(idx2));
}

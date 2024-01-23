/***************************************************************/
//
//Copyright (c) 2022 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/
#include "results.h"
#include <iomanip>
#include <sstream>
#include <fstream>
using namespace opal;
using namespace optix;


void ResultRecord::insertHitRecord(HitRecord h) {
	hitRecords.push_back(h);
}
void ResultRecord::clearHits() {
	hitRecords.clear();
}
std::string HitRecord::toString() const  {
	std::stringstream os;
        os<<E.x<<","<<E.y<<","<<Ex.x<<","<<Ex.y<<","<<Ey.x<<","<<Ey.y<<","<<Ez.x<<","<<Ez.y<<","<<unfoldedPath<<","<<directionOfArrival.x<<","<<directionOfArrival.y<<","<<directionOfArrival.z<<","<<directionOfDeparture.x<<","<<directionOfDeparture.y<<","<<directionOfDeparture.z<<","<<diffraction;
	return os.str();
    
}
void HitCollector::insertHit(int txId, int rxId, HitRecord h) {
	auto tx=hits.find(txId);
	if (tx==hits.end()) {
		//New transmitter
		HitEntry* he = new HitEntry();
		std::vector<HitRecord>* v=new std::vector<HitRecord>();
		//std::cout<<"Insert new tx "<<txId<<"to rx "<<" "<<rxId<<h.toString()<<std::endl;
		v->push_back(h);
		he->insert(std::pair<int,std::vector<HitRecord>*>(rxId,v));
		hits.insert(std::pair<int,HitEntry*>(txId,he));

	} else {
		//tx exists
		auto rx=tx->second->find(rxId);
		if (rx!=tx->second->end()) {
			std::vector<HitRecord>* v=(rx->second);
			//Receiver exists, add record
			v->push_back(h);
		} else {
			std::vector<HitRecord>* v=new std::vector<HitRecord>();
			v->push_back(h);
			tx->second->insert(std::pair<int,std::vector<HitRecord>*>(rxId,v));
			

		}

	}
}
std::vector<HitRecord> HitCollector::getRecords(int txId, int rxId) {
	if (txId==rxId) {
                //We cannot have self hits, return empty vecotr
                std::vector<HitRecord> empty;
                return empty;
        }
	auto tx=hits.find(txId);
	if (tx!=hits.end()) {
		auto rx=tx->second->find(rxId);
		if (rx!=tx->second->end()) {
			//auto v=(*(rx->second));
			//for (int i =0; i< v.size(); ++i) {
			//	std::cout<<"Records for  tx "<<txId<<"to rx "<<rxId<<" "<<v[i].toString()<<std::endl;
			//	
			//}
			return (*(rx->second));
		} else {
			std::cout << "HitCollector::getRecords(): Not found receiver " << rxId<< "for transmitter" <<txId<< std::endl;
			throw  opal::Exception("HitCollector::getRecords(): rx not found");
		}
	} else {
		std::cout << "HitCollector::getRecords(): Not found transmitter " << txId<< std::endl;
		throw  opal::Exception("HitCollector::getRecords(): tx not found");
	}

}
HitCollector::~HitCollector() {
	clear();
}
void HitCollector::clear() {
	for (auto const& he : hits) {
		for (auto const& h: (*he.second)) {
			delete h.second;
		}
		delete he.second;
	}
	hits.clear();
}
void ResultRecord::setField(optix::float2 E, optix::float2 Ex,optix::float2 Ey,optix::float2 Ez,unsigned int index, int rxId, int txId, optix::float3 position, float radius, optix::float3 origin, float eA, float txPower, int refHits, int difHits, float frequency ) 
{
	this->E=E;
	this->Ex=Ex;
	this->Ey=Ey;
	this->Ez=Ez;
	this->index=index;
	this->rxId=rxId; 
	this->txId=txId;
	this->position=position;
	this->radius=radius;
	this->origin=origin; 
	this->eA=eA; 
	this->txPower=txPower; 
	this->refHits=refHits;
	this->difHits=difHits;
	this->d=length(origin-position);
	this->frequency=frequency;
	this->power=eA*((E.x*E.x) + (E.y*E.y))*txPower;

}
std::string ResultRecord::toString() const  {
	std::stringstream os;
		os<<std::setprecision(10)<<"TPR\t"<<txId<<"\t"<<rxId<<"\t"<<power<<"\t"<<position.x<<"\t"<<position.y <<"\t"<<position.z<<"\t"<<refHits<<"\t"<<difHits<<"\t"<<radius<<"\t"<<d<<"\t"<<origin.x<<"\t"<<origin.y<<"\t"<<origin.z<<"\t"<<eA <<"\t"<<txPower<<"\t"<<frequency<< std::endl;
		os<<std::setprecision(10)<<"TER\t"<<txId<<"\t"<<rxId<<"\t"<<E.x<<"\t"<<E.y<<"\t"<<position.x<<"\t"<<position.y <<"\t"<<position.z<<"\t"<<refHits<<"\t"<<difHits<<"\t"<<radius<<"\t"<<d<<"\t"<<origin.x<<"\t"<<origin.y<<"\t"<<origin.z<<"\t"<<eA<<"\t"<<txPower<<"\t"<<frequency<< std::endl;
		os<<"TERX\t"<<txId<<"\t"<<rxId<<"\t"<<Ex.x<<"\t"<<Ex.y<<"\t"<<position.x<<"\t"<<position.y<<"\t"<<position.z<<"\t"<<refHits<<"\t"<<difHits<<"\t"<<radius<<"\t"<<d<<"\t"<<origin.x<<"\t"<<origin.y<<"\t"<<origin.z<<"\t"<<txPower<<"\t"<<frequency<< std::endl;
		os<<"TPRX\t"<<txId<<"\t"<<rxId<<"\t"<<(10*log10(dot(Ex,Ex)))<<"\t"<<dot(Ex,Ex)<<"\t"<<position.x<<"\t"<<position.y<<"\t"<<position.z<<"\t"<<refHits<<"\t"<<difHits<<"\t"<<radius<<"\t"<<d<<"\t"<<origin.x<<"\t"<<origin.y<<"\t"<<origin.z<< "\t"<<txPower<<"\t"<<frequency<< std::endl;

		os<<"TERY\t"<<txId<<"\t"<<rxId<<"\t"<<Ey.x<<"\t"<<Ey.y<<"\t"<<position.x<<"\t"<<position.y<<"\t"<<position.z<<"\t"<<refHits<<"\t"<<difHits<<"\t"<<radius<<"\t"<<d<<"\t"<<origin.x<<"\t"<<origin.y<<"\t"<<origin.z<<"\t"<<txPower<< "\t"<<frequency<< std::endl;

		os<<"TPRY\t"<<txId<<"\t"<<rxId<<"\t"<<(10*log10(dot(Ey,Ey)))<<"\t"<<dot(Ey,Ey)<<"\t"<<position.x<<"\t"<<position.y<<"\t"<<position.z<<"\t"<<refHits<<"\t"<<difHits<<"\t"<<radius<<"\t"<<d<<"\t"<<origin.x<<"\t"<<origin.y<<"\t"<<origin.z<<"\t"<<txPower<<"\t"<<frequency<<  std::endl;

		os<<"TERZ\t"<<txId<<"\t"<<rxId<<"\t"<<Ez.x<<"\t"<<Ez.y<<"\t"<<position.x<<"\t"<<position.y<<"\t"<<position.z<<"\t"<<refHits<<"\t"<<difHits<<"\t"<<radius<<"\t"<<d<<"\t"<<origin.x<<"\t"<<origin.y<<"\t"<<origin.z<<"\t"<<txPower<< "\t"<<frequency<< std::endl;

		os<<"TPRZ\t"<<txId<<"\t"<<rxId<<"\t"<<(10*log10(dot(Ez,Ez)))<<"\t"<<dot(Ez,Ez)<<"\t"<<position.x<<"\t"<<position.y<<"\t"<<position.z<<"\t"<<refHits<<"\t"<<difHits<<"\t"<<radius<<"\t"<<d<<"\t"<<origin.x<<"\t"<<origin.y<<"\t"<<origin.z<<"\t"<<txPower<<"\t"<<frequency<< std::endl;


	
	return os.str();
}
std::ostream& operator<<(std::ostream& os, const opal::ResultRecord& r)  {
	os<< r.toString();
	return os;
}

ResultReport::ResultReport() {
}
void ResultReport::insert(ResultRecord r) {
	results.push_back(r);
}
std::string ResultReport::toString() const {
	std::stringstream os;
	for (auto r : results) {
		os<<r.toString();
	}
	return os.str();
}
void ResultReport::merge(ResultReport& other) {
	auto o=other.getVector();
	results.insert(results.end(),o.begin(), o.end());
}
void ResultReport::sort() {
	std::sort(results.begin(), results.end());
}	
void ResultReport::sortByRx() {
	byRx s;
	std::sort(results.begin(), results.end(), s );
}	
std::ostream& operator<<(std::ostream& os, const opal::ResultReport& rr) {
	os<<rr.toString();
	return os;
}

void ResultReport::toCSV(std::string file) {
	std::ofstream myfile;
	myfile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
	try {
		if (results.size()>0) {
			myfile.open(file);
			myfile << "txId,rxId,V.real,V.im,power,Ex.real,Ex.im,Px,Ey.real,Ey.im,Py,Ez.real,Ez.im,Pz,position.x,position.y,position.z,refHits,difHits,radius,distance,txpos.x,txpos.y,txpos.z,txPower,frequency" << std::endl;
			for (auto r : results) {
				myfile << std::setprecision(10)<<  r.txId << "," << r.rxId << "," <<r.E.x<<","<<r.E.y<<","<<r.power<<","<< r.Ex.x << "," << r.Ex.y << "," << dot(r.Ex, r.Ex) << "," << r.Ey.x << "," << r.Ey.y << "," << dot(r.Ey, r.Ey) << "," << r.Ez.x << "," << r.Ez.y << "," << dot(r.Ez, r.Ez) << "," << r.position.x << "," << r.position.y << "," << r.position.z << "," << r.refHits << "," << r.difHits << "," << r.radius << "," << r.d << "," << r.origin.x << "," << r.origin.y << "," << r.origin.z << "," << r.txPower << "," <<r.frequency<<std::endl;

			}
			//myfile << "E.real,E.im,power,refHits,difHits" << std::endl;

			//for (auto r : results) {
			//	myfile << std::setprecision(10)  << r.E.x << "," << r.E.y << "," << r.power << "," << r.refHits << "," << r.difHits << std::endl;
			//}


			myfile.close();
		} else {
			std::cout<<"WARNING: no results recorded. Not writing CSV results file"<<std::endl;
		}
	}
	catch (std::ifstream::failure e) {
		std::cerr << "ResultReport::toCSV(): Exception opening/reading/closing file "<<file<<" \n"; 

	}
}


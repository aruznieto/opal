/***************************************************************/
//
//Copyright (c) 2019 Esteban Egea-Lopez http://girtel.upct.es/~simulations
//
/**************************************************************/


#ifndef RESULTS_H
#define RESULTS_H
#include "Opal.h"
#include <vector>
#include <iostream>
#include <tuple>
#include <algorithm>
namespace opal {



	//Record information about individual hits for on one transmission
	class HitRecord {
		public:
			bool diffraction;
			optix::float2 E;
			optix::float2 Ex;
			optix::float2 Ey;
			optix::float2 Ez;
			optix::float3 directionOfArrival;
			optix::float3 directionOfDeparture;
			float unfoldedPath;
			//Aditional info can be added here
			std::string toString() const ;

	};
	class HitCollector {
		protected:
			typedef std::map<int,std::vector<HitRecord>*> HitEntry; //Set of hits for receiver
			typedef std::map<int, HitEntry*> HitMap;
			HitMap hits;
		public:
			virtual ~HitCollector();
			void insertHit(int txId, int rxId, HitRecord h);
			std::vector<HitRecord>  getRecords(int txId, int rxId) ;
			void clear() ;
	};
	class ResultRecord  {
		public:
			unsigned int index;
			int rxId;
			int txId;
			optix::float3 position;
			float radius;
			optix::float3 origin;
			float eA;
			float txPower;
			int refHits;
			int difHits;
			optix::float2 E;
			optix::float2 Ex;
			optix::float2 Ey;
			optix::float2 Ez;
			float frequency;	
			float power;
			float d;
		 	void setField(optix::float2 E, optix::float2 Ex, optix::float2 Ey, optix::float2 Ez, unsigned int index, int rxId, int txId, optix::float3 positio, float radius, optix::float3 origin, float eA, float txPower, int refHits, int difHits, float frequency );	
			//Default sort, first by txId and then by rxId
			inline  bool operator<(const ResultRecord &h) const {
				if (txId==h.txId) {
					return (rxId<h.rxId);
				} else {
					return (txId<h.txId);
				}	
			};
			std::string toString() const ;
		        std::vector<HitRecord> hitRecords;
		        void insertHitRecord(HitRecord h) ;	
			void clearHits();
			//friend std::ostream& operator<<(std::ostream& os, const ResultRecord& r) ;
		
	};
	//Sort options
	//sort by id of receiver
	struct byRx {
			bool operator()(ResultRecord a, ResultRecord b) const { 
			return std::tie(a.rxId,a.txId) < std::tie(b.rxId,b.txId) ; }
	} ;
	//sort by z position of receiver and then receiver id
	struct byZRx {
			bool operator()(ResultRecord a, ResultRecord b) const { 
			return std::tie(a.position.z,a.rxId,a.txId) < std::tie(b.position.z,b.rxId,b.txId) ; }
	} ;

	//sort by id of receiver and then receiver z position
	struct byRxZ {
			bool operator()(ResultRecord a, ResultRecord b) const { 
			return std::tie(a.rxId,a.position.z,a.txId) < std::tie(b.rxId,b.position.z,b.txId) ; }
	} ;
	class ResultReport {
		//This may be done as templates
		public:
			std::vector<ResultRecord> results;
			ResultReport() ;
			void insert(ResultRecord r);
			void clear() {results.clear();} ;
			unsigned int size() const { return results.size();};
			std::vector<ResultRecord>& getVector() { return results;};	
			std::string toString() const ;
 			void toCSV(std::string file);
			void sort();
			void sortByRx();
			void merge(ResultReport& other);
			template <typename T> 
			void sortBy(T t) { 
				std::sort(results.begin(), results.end(),  t);
			};
	};
}
#endif


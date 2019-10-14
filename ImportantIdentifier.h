#ifndef IMPORTANT_IDENTIFIER_H
#define IMPORTANT_IDENTIFIER_H

#include "GiscupParser.h"
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <cstring>
#include <cassert>
#include <iostream>
#include "utils.h"



template <class GlobalIdType>
class ImportantIdentifier;

/*
Here we define three ways to store the global ids:

We also define the ImportantIdentifier class, that represents important tokens (vidaId, fromId, rows, controllers, etc) found in the input JSON file
*/


/*
-- GlobalIdTypeCharPtrs: the ids are represented by two char pointers pointing to the begin and end of each global id. The disadvantage
is that the entire json file must be kept in memory during all the lifetime of the program (the pointers point to the regions
of the file containing the global ids). This representation cannot be used depending on the kind of the parser we employ (one of the parsers
read the file in small blocks and, thus, this kind of representation cannot be used with it).
*/
class GlobalIdTypeCharPtrs {
public:
	GlobalIdTypeCharPtrs(const char *idStartPos_=NULL,const char* idEndPos_=NULL):idStartPos(idStartPos_),idEndPos(idEndPos_) {}
		

	size_t stringIdSize() const {return idEndPos-idStartPos+1;}

	bool compareLess(const GlobalIdTypeCharPtrs &b) const {
		const GlobalIdTypeCharPtrs &a = *this;

		int as= a.stringIdSize();
  	int bs= b.stringIdSize();
  	if(as!=bs) { 
  		if(as<bs) {
  			int diff = strncmp(a.idStartPos,b.idStartPos,a.stringIdSize());
  			if(diff!=0) return diff<0;
  			else return true; //a is smaller than b... if the common prefix is equal --> a is lexicographicallty smaller
  		} else {
  			int diff = strncmp(a.idStartPos,b.idStartPos,b.stringIdSize());
  			if(diff!=0) return diff<0;
  			else return false; //a is greater than b... if the common prefix is equal --> a is lexicographicallty greater
  		}  		
  	} 
  	assert(a.stringIdSize()==b.stringIdSize());
  	return strncmp(a.idStartPos,b.idStartPos,a.stringIdSize())<0;
	}

	bool compareEq(const GlobalIdTypeCharPtrs &b) const {
		const GlobalIdTypeCharPtrs &a = *this;

		int as= a.stringIdSize();
  	int bs= b.stringIdSize();
  	return (as==bs) && strncmp(a.idStartPos,b.idStartPos,a.stringIdSize())==0;
	}

	void print(std::ostream &out) const {
		out.write(idStartPos,idEndPos-idStartPos+1);	
	}

	std::string toString() const {
		return std::string(idStartPos,idEndPos-idStartPos+1);
	}
private:
	const char *idStartPos,*idEndPos;
};


/*
This is the simplest representation:
-- GlobalIdTypeString : the global ids are represented by strings (e.g.: "{3AADAD81-B0C6-83CF-611B-C7CA366D9B4A}"). This has some performance
overheads because we have to store a string with several bytes (for each global id) and, also, because this leads to a lot of
small allocations in the heap.
*/
class GlobalIdTypeString {
public:
	//pointers to the first and last char of the id
	GlobalIdTypeString(const char *idStartPos=NULL,const char* idEndPos=NULL) {
		if(idStartPos)
			guid = string(idStartPos,idEndPos+1);
	}
		

	size_t stringIdSize() const {return guid.size();}

	bool compareLess(const GlobalIdTypeString &b) const {
		return guid < b.guid;
	}

	bool compareEq(const GlobalIdTypeString &b) const {
		return guid==b.guid;
	}

	void print(std::ostream &out) const {
		out.write(&(guid[0]),stringIdSize());	
	}

	std::string toString() const {
		return guid;
	}
private:
	string guid;
};




//Converts an hexadecimal (integer) digit into a char
//0 -> 0
//1 -> 1
//10 -> A
//15 -> F
char toChar(int i) {
	if(i<10) return i+'0';
	return i-10+'A';
}

//Converts an hexadecimal digit (char) into an integer (from 0 to 15)
inline int decodeHexChar(const char c) {
	return 9*(c>>6) + (c&0xf); //bitwise trick to convert an hex char into an integer...
}



//decodes a string (GUID) represented as {3AADAD81-B0C6-83CF-611B-C7CA366D9B4A}, converting it into two 64-bit unsigned long longs
static const int startVA[] = {0,4,9,14}; //start position of each 4 characters (16 bits) of the GUID string (ignoring the { and }) -- (here we consider the first 64 bits of the string)
static const int startVB[] = {19,24,28,32}; //start position of each 4 characters (16 bits) of the GUID string (ignoring the { and }) -- (here we consider the last 64 bits of the string)
pair<unsigned long long,unsigned long long> decode(const char *str, const char *end) {
		if(end-str!=35) throw -1;
		unsigned long long partAv[4] = {0}; //first part of the ID (first 64 bits)
		//For performance, we created a code that is easier to vectorize... (but a little less easy to read)
		//we decode 4 characters each time...

		//decodes first half
		for(int j=0;j<4;j++)  {
			for(int k=0;k<4;k++) 
				partAv[k] <<= 4;
			for(int k=0;k<4;k++) 
				partAv[k] |= decodeHexChar(str[j + startVA[k]]);			
		}

		//second half...
		unsigned long long partBv[4] = {0}; //second part of the ID (last 64 bits)
		for(int j=0;j<4;j++)  {
			for(int k=0;k<4;k++) 
				partBv[k] <<= 4;
			for(int k=0;k<4;k++) 
				partBv[k] |= decodeHexChar(str[j + startVB[k]]);
		}


		return make_pair((partAv[0]<<48|partAv[1]<<32|partAv[2]<<16|partAv[3]) ,  (partBv[0]<<48|partBv[1]<<32|partBv[2]<<16|partBv[3]));
}


/*
-- GlobalIdTypeGUID: assumming the global ids always have the format shown above ("{3AADAD81-B0C6-83CF-611B-C7CA366D9B4A}"), we can
codify it using two 64-bit integers (each number/letter is an hexadecimal digit). 
*/
class GlobalIdTypeGUID {
public:
	GlobalIdTypeGUID(const char *idStartPos=NULL,const char* idEndPos=NULL) {
		if(!idStartPos) return;
		guid = decode(idStartPos+1,idEndPos-1); //ignore the { and }
	}
		

	size_t stringIdSize() const {return 38;}

	bool compareLess(const GlobalIdTypeGUID &b) const {
		return guid < b.guid;
	}

	bool compareEq(const GlobalIdTypeGUID &b) const {
		return guid==b.guid;
	}

	//prints the identifier using the format "{06A51F9A-D129-C762-CDAE-209E54B4603B}"
	void print(std::ostream &out) const {
		//06A51F9A-D129-C762-CDAE-209E54B4603B
		
		char temp[38];
		temp[0] =  '{';
		temp[9] = '-';
		temp[14] = '-';
		temp[19] = '-';
		temp[24] = '-';
		temp[37] = '}';

		unsigned long long f = guid.first;
		int i=18;
		for(;i!=14;i--,f>>=4) temp[i] = toChar(f&15); 
		i--;
		for(;i!=9;i--,f>>=4) temp[i] = toChar(f&15); 
		i--;
		for(;i>0;i--,f>>=4) temp[i] = toChar(f&15); 
		i--;

		f = guid.second;
		i=36;
		for(;i!=24;i--,f>>=4) temp[i] = toChar(f&15); 
		i--;
		for(;i!=19;i--,f>>=4) temp[i] = toChar(f&15); 
		i--;
		out.write(temp,38);	
		return;
	}


private:
	friend ImportantIdentifier<GlobalIdTypeGUID>;
	pair<unsigned long long, unsigned long long> guid;
};



template <class GlobalIdType>
struct ImportantIdentifierHasher;

template <class GlobalIdType>
struct ImportantIdentifierEqComparisonString;

template <class T>
struct ImportantIdentifierLessThanComparisonString;

template <class GlobalIdType>
struct ImportantIdentifierHashLessThanComparison; 










/*
Represent an important token found in the input file.

For example, if we are using the GlobalIdTypeGUID type of identification, a fromGlobalId token will:
-Have a GlobalIdTypeGUID object representing the from global id.
-Have the type:  FROM_IDENTIFIER
-Have a hash: we pre-compute the hash of the global ids so that we can efficiently add the tokens to hash tables, get unique tokens (for example, when 
writing the output), etc.
*/

template <class GlobalIdType>
class ImportantIdentifier {
public:
	typedef ImportantIdentifierHasher<GlobalIdType>  ImportantIdentifierHasherType;
	typedef ImportantIdentifierEqComparisonString<GlobalIdType>  ImportantIdentifierEqComparisonStringType;
	typedef ImportantIdentifierLessThanComparisonString<GlobalIdType>  ImportantIdentifierLessThanComparisonStringType;
	typedef ImportantIdentifierHashLessThanComparison<GlobalIdType>  ImportantIdentifierHashLessThanComparisonType;

	ImportantIdentifier(int type_): id(-1){
		setType(type_);
	}

	ImportantIdentifier(int type_, const char *idStartPos_,const char* idEndPos_);

	ImportantIdentifier(): globalId(NULL,NULL), id(-1){
		setType(-1);
	}
	void setType(int type_) {
		type = type_;
	}
	void print(std::ostream &out) const {
		globalId.print(out);
	}
	std::string toString() const {
		return globalId.toString();
	}

	

	//vertices are labeled 0...nv-1 (nv is the number of vertices)
	//edges are labeled 0...ne-1 (ne is the number of edges)
	int id; //internal id of the identifier 

	size_t stringIdSize() const {return globalId.stringIdSize();}

	//type of token
	const static int ROWS_IDENTIFIER = 0;
	const static int CONTROLLERS_IDENTIFIER = 1;
	const static int GLOBAL_ID_IDENTIFIER = 2;
	const static int FROM_IDENTIFIER = 3;
	const static int VIA_IDENTIFIER = 4;
	const static int TO_IDENTIFIER = 5;
	const static int VERTEX_STARTING_POINT_IDENTIFIER = 6;
	const static int EDGE_STARTING_POINT_IDENTIFIER = 7;
	int type;

	//For performance, we also store the hash of the tokens
	unsigned long long hash;	
private:	
	GlobalIdType globalId;	

	friend ImportantIdentifierHasherType;
	friend ImportantIdentifierEqComparisonStringType;
	friend ImportantIdentifierHashLessThanComparisonType;
	friend ImportantIdentifierLessThanComparisonStringType;
};

//When we use the type GlobalIdTypeGUID, the hash is computed based on the 2 unsigned long long representing the number
//(this is faster because it avoids processing all the chars of the ID again...s)
template <>
inline ImportantIdentifier<GlobalIdTypeGUID>::ImportantIdentifier(int type_, const char *idStartPos_,const char* idEndPos_): globalId(idStartPos_,idEndPos_),id(-1){
	setType(type_);
	hash = globalId.guid.first + globalId.guid.second;
}

template <class GlobalIdType>
inline ImportantIdentifier<GlobalIdType>::ImportantIdentifier(int type_, const char *idStartPos_,const char* idEndPos_): globalId(idStartPos_,idEndPos_),id(-1){
	setType(type_);

	hash = 5381;
	for(const char *ptr=idStartPos_;ptr<=idEndPos_;ptr++) {
		const char c = *ptr;
		hash = ((hash << 5) + hash) + ((int)c);
	}
}

//Comparison functors employed to compare identifiers

//Gets the hash of an identifier
template <class GlobalIdType>
struct ImportantIdentifierHasher {
  std::size_t operator()(const ImportantIdentifier<GlobalIdType>& identifier) const {
    return identifier.hash;
  }
};

//Compares the global id for equality 
template <class GlobalIdType>
struct ImportantIdentifierEqComparisonString {
  bool operator()(const ImportantIdentifier<GlobalIdType>& a,const ImportantIdentifier<GlobalIdType>& b) const {	  
  	if(a.hash!=b.hash) return false;
  	
    return a.globalId.compareEq(b.globalId);
  }
};

//Compares the global id for <= 
template <class GlobalIdType>
struct ImportantIdentifierLessThanComparisonString {
  std::size_t operator()(const ImportantIdentifier<GlobalIdType>& a,const ImportantIdentifier<GlobalIdType>& b) const {
  	return a.globalId.compareLess(b.globalId);
  }
	};

//Compares the hash of the global ids for <=
//We sometimes need this...
template <class GlobalIdType>
struct ImportantIdentifierHashLessThanComparison {
	bool operator()(const ImportantIdentifier<GlobalIdType>& a,const ImportantIdentifier<GlobalIdType>& b) const {
			unsigned long long ahash = a.hash;
			unsigned long long bhash = b.hash;
			return ahash < bhash;
	}
};






#endif
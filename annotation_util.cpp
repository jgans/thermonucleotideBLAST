#include "annotation.h"

#include <sstream>
#include <string.h>

using namespace std;

#define	MAX_LINE_LEN	1024

unsigned long int line_number = 1;

bool is_eol(char m_c);

// Note that genbank entries are 1's based (not zero based like asn.1
// entires).
bool read_range(gzFile m_fin, pair<unsigned int, unsigned int> &m_range, 
				list< pair<unsigned int, unsigned int> > &m_seg_list)
{
	unsigned int i = 0;
	char buffer[MAX_LINE_LEN];

	if( (gzgets(m_fin, buffer, MAX_LINE_LEN) == NULL) || !strip_eol(buffer, MAX_LINE_LEN) ){
		throw __FILE__ ":read_range: Unable to read line";
	}

	++line_number;
	
	unsigned int len = strlen(buffer);

	// Make the range parsing robust to Mac/PC/Unix EOL differences
	if((len > 0) && (buffer[len - 1] == '\r') ){
		len--;
	}
	
	// Skip to the start of the range entry.
	while(isspace(buffer[i])){
		i++;
	}
	
	// For now, ignore "<" and ">"
	if((buffer[i] == '<') || (buffer[i] == '>')){
		i++;
	}

	// Is this a basic range (i.e. xxx..yyy) ?
	if( isdigit(buffer[i]) ){

		// Read a basic range: The first entry ...
		list<char> number;

		while( isdigit(buffer[i]) ){
		
			number.push_back(buffer[i] - '0');
			i++;
		}
	
		// Convert this list to an zero based base location
		m_range.first = list_to_int(number) - 1;

		if(i == len){

			// We have found a single number entry for the
			// range (i.e. NC_005816.gbk)
			m_range.second = m_range.first;

			// Not a complement
			return false;
		}

		// The second entry. Read the spacers ".."
		while((i != len) && !isdigit(buffer[i])){
			i++;
		}

		while((i != len) && isdigit(buffer[i])){
			number.push_back(buffer[i] - '0');
			i++;
		}
		
		// Convert this list to an zero based base location
		m_range.second = list_to_int(number) - 1;

		// Not a complement
		return false;
	}

	// Is this range a simple complement?
	if(buffer[i] == 'c'){
	
		unsigned int j = i + 11; // + strlen( "complement(" )
		
		// For now, ignore "<" and ">"
		if((buffer[j] == '<') || (buffer[j] == '>')){
			j++;
		}

		// Look -- not checking the bounds on buffer for speed reasons
		if( isdigit(buffer[j]) ){
		
			// Read a basic range: The first entry ...
			list<char> number;

			while(isdigit(buffer[j])){
				number.push_back(buffer[j] - '0');
				j++;
			}
			
			// Convert this list to an zero based base location
			m_range.first = list_to_int(number) - 1;

			if(j == len - 1 /* skip the closing ')' of complement(...) */){

				// We have found a single number entry for the
				// range (i.e. NC_005816.gbk)
				m_range.second = m_range.first;

				// Is a complement
				return true;
			}

			// The second entry. Read the spacers ".."
			while( !isdigit(buffer[j]) ){
				j++;
			}

			while( isdigit(buffer[j]) ){
			
				number.push_back(buffer[j] - '0');
				j++;
			}
			
			// Convert this list to an zero based base location
			m_range.second = list_to_int(number) - 1;

			// Is a complement
			return true;
		}
	}
	
	// We are now into the complicated ranges! i.e. join and complement(join( .. ) )
	// First, make sure that we have the entire range!
	unsigned int left_paren = 0;
	unsigned int right_paren = 0;
	unsigned int j;

	for(j = i;buffer[j] != '\0';j++){
		if(buffer[j] == '('){
			left_paren ++;
		}

		if(buffer[j] == ')'){
			right_paren ++;
		}
	}

	while(left_paren != right_paren){
		
		char tmp[MAX_LINE_LEN];

		if( (gzgets(m_fin, tmp, MAX_LINE_LEN) == NULL) || !strip_eol(tmp, MAX_LINE_LEN) ){
			throw __FILE__ ":read_range: Unable to read tmp line";
		}

		line_number ++;

		if( ( strlen(buffer) + strlen(tmp) ) >= MAX_LINE_LEN){
			throw __FILE__ ":read_range: Ran out of room to store range!";
		}

		strcat(buffer, tmp);
		
		left_paren = right_paren = 0;

		for(unsigned int j = i;buffer[j] != '\0';j++){

			if(buffer[j] == '('){
				left_paren ++;
			}

			if(buffer[j] == ')'){
				right_paren ++;
			}
		}
	}

	bool is_complement = false;
	//bool is_join = false;

	// Now test for complement
	if(buffer[i] == 'c'){
		is_complement = true;

		i += 11; // + strlen( "complement(" )
	}

	// Read the join info
	if(buffer[i] == 'j'){
		//is_join = true;

		i += 5; // + strlen( "join(" )
	}
	else{
		if(buffer[i] == 'o'){
			//is_join = true;

			i += 6; // + strlen( "order(" )
		}
		else{
			throw error_msg(":read_range: Did not find expected \"join\"");
		}
	}
	
	// Test for "complement" *again* to catch the extreamly rare "join(complement" -- see 
	// NC_005213.gbk for an example.
	if(buffer[i] == 'c'){
		is_complement = true;

		i += 11; // + strlen( "complement(" )
	}
	
	// For now, ignore "<" and ">"
	if((buffer[i] == '<') || (buffer[i] == '>')){
		i++;
	}

	j = strlen(buffer);

	while(i < j){
		// Read a range: The first entry ...
		list<char> number;
		
		pair<unsigned int, unsigned int> tmp_range;

		while((i < j) && isdigit(buffer[i])){
			number.push_back(buffer[i] - '0');
			i++;
		}
		
		// Convert this list to an zero based base location
		tmp_range.first = list_to_int(number) - 1;
		
		// Is there a second value in this range, or is it a single
		// base entry?
		bool is_single_base = false;

		while((i < j) && !isdigit(buffer[i])){
			if((buffer[i] == ',') || (buffer[i] == ')')){
				is_single_base = true;
			}

			i++;
		}

		// Is there a second entry to read?
		if(is_single_base){
			tmp_range.second = tmp_range.first;
		}
		else{
			// The second entry
			while((i < j) && isdigit(buffer[i])){
				number.push_back(buffer[i] - '0');
				i++;
			}
			
			// Convert this list to an zero based base location
			tmp_range.second = list_to_int(number) - 1;
		}

		m_seg_list.push_back(tmp_range);

		// Skip all non-digit characters
		while((i < j) && !isdigit(buffer[i])){
			i++;
		}
	}

	if(m_seg_list.empty()){
		throw error_msg("read_range: Unable to read join(...)");
	}

	// Sort the seg list (important for annotations that span the origin)
	m_seg_list.sort();

	m_range.first = m_seg_list.front().first;
	m_range.second = m_seg_list.back().second;

	return is_complement;

	
}

// Convert the input list of char's into an unsinged integer.
// The list is exhausted in the process
unsigned int list_to_int(list<char> &m_number)
{
#ifdef _DEBUG
	if(m_number.empty()){
		throw error_msg("Empty number list!");
	}
#endif // _DEBUG

	unsigned int num = m_number.back();
			
	m_number.pop_back();

	int power = 10;

	while(m_number.empty() == false){
		num += power*m_number.back();

		power *= 10;
		m_number.pop_back();
	}
			
	return num;
}

char *error_msg(const char *m_error)
{
	static char error[1024];
	
	stringstream ssout;
	
	ssout << m_error << " : line #" << line_number;
	
	strcpy( error, ssout.str().c_str() );
	
	return error;
}

bool is_eol(char m_c)
{
	return (m_c == '\n') || (m_c == '\r');
}

// Return true if we successfully found and removed any trailing eol
// characters (i.e., '\n' or '\r')
bool strip_eol(char *m_buffer, const size_t &m_max_len)
{
	size_t len = strnlen(m_buffer, m_max_len);

	if(len == m_max_len){
		return false;
	}

	bool found_eol = false;

	while( (len > 0) && is_eol(m_buffer[len - 1]) ){

		found_eol = true;
		m_buffer[len - 1] = '\0';
		--len;
	}

	return found_eol;
}
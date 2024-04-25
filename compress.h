#ifndef __COMPRESS_DNA
#define __COMPRESS_DNA

#include <string>

// Compress DNA sequences using the restricted alphabet defined in seq.h and "hand carved" Huffman encoding
// (based on estimated symbol frequencies, not actual measurements).

// Real bases (these are most common, 3 bits per base)
#define     HUFFMAN_A       0b000
#define     HUFFMAN_T       0b001
#define     HUFFMAN_G       0b010
#define     HUFFMAN_C       0b011
#define     HUFFMAN_END     0b100       // Every compressed sequence requries a termination symbol

// Less common symbols, 6 bit
#define     HUFFMAN_N       0b101000       // The completely degenerate base is often used to mask unknown sequences 
#define     HUFFMAN_3       0b101001
#define     HUFFMAN_5       0b101010
#define     HUFFMAN_SPACE   0b101011
#define     HUFFMAN_EOL     0b101100
#define     HUFFMAN_PRIME   0b101101
#define     HUFFMAN_PIPE    0b101110
#define     HUFFMAN_COLON   0b101111

// Rarely used degenerate bases and other symbols, 7 bits per symbol
#define     HUFFMAN_GAP     0b1100000
#define     HUFFMAN_M       0b1100001
#define     HUFFMAN_R       0b1100010
#define     HUFFMAN_S       0b1100011
#define     HUFFMAN_V       0b1100100
#define     HUFFMAN_W       0b1100101
#define     HUFFMAN_Y       0b1100110
#define     HUFFMAN_H       0b1100111
#define     HUFFMAN_K       0b1101000
#define     HUFFMAN_D       0b1101001
#define     HUFFMAN_B       0b1101010
#define     HUFFMAN_I       0b1101011
// Lower case bases are used to indicate binding sites, also 7 bits per symbol
#define     HUFFMAN_a       0b1101100
#define     HUFFMAN_t       0b1101101
#define     HUFFMAN_g       0b1101110
#define     HUFFMAN_c       0b1101111
#define     HUFFMAN_n       0b1110000
#define     HUFFMAN_m       0b1110001
#define     HUFFMAN_r       0b1110010
#define     HUFFMAN_s       0b1110011
#define     HUFFMAN_v       0b1110100
#define     HUFFMAN_w       0b1110101
#define     HUFFMAN_y       0b1110110
#define     HUFFMAN_h       0b1110111
#define     HUFFMAN_k       0b1111000
#define     HUFFMAN_d       0b1111001
#define     HUFFMAN_b       0b1111010
#define     HUFFMAN_i       0b1111011

// Compress a DNA sequence
inline std::string deflate_dna_seq(const std::string &m_seq)
{
    std::string ret;

    char c = 0x0;
    uint8_t bit = 0;

    #define PACK(VALUE, NUM_BITS) \
        for(int8_t index = 0;index < NUM_BITS;++index){ \
            c = (c << 1) | ( 1 & (VALUE >> (NUM_BITS - 1 - index) ) ); \
            ++bit; \
            if(bit == 8){\
                ret.push_back(c); \
                c = 0x0; \
                bit = 0; \
            } \
        }

    for(std::string::const_iterator i = m_seq.begin();i != m_seq.end();++i){

        switch(*i){
            case 'A':
                PACK(HUFFMAN_A, 3);
                break;
            case 'T':
                PACK(HUFFMAN_T, 3);
                break;
            case 'G':
                PACK(HUFFMAN_G, 3);
                break;
            case 'C':
                PACK(HUFFMAN_C, 3);
                break;
            case 'N':
                PACK(HUFFMAN_N, 6);
                break;
             case '3':
                PACK(HUFFMAN_3, 6);
                break;
            case '5':
                PACK(HUFFMAN_5, 6);
                break;
            case ' ':
                PACK(HUFFMAN_SPACE, 6);
                break;
            case '\n':
                PACK(HUFFMAN_EOL, 6);
                break;
            case '\'':
                PACK(HUFFMAN_PRIME, 6);
                break;
            case '|':
                PACK(HUFFMAN_PIPE, 6);
                break;
            case ':':
                PACK(HUFFMAN_COLON, 6);
                break;
            case '-':
                PACK(HUFFMAN_GAP, 7);
                break;
            case 'M':
                PACK(HUFFMAN_M, 7);
                break;
            case 'R':
                PACK(HUFFMAN_R, 7);
                break;
            case 'S':
                PACK(HUFFMAN_S, 7);
                break;
            case 'V':
                PACK(HUFFMAN_V, 7);
                break;
            case 'W':
                PACK(HUFFMAN_W, 7);
                break;
            case 'Y':
                PACK(HUFFMAN_Y, 7);
                break;
            case 'H':
                PACK(HUFFMAN_H, 7);
                break;
            case 'K':
                PACK(HUFFMAN_K, 7);
                break;
            case 'D':
                PACK(HUFFMAN_D, 7);
                break;
            case 'B':
                PACK(HUFFMAN_B, 7);
                break;
            case 'I':
                PACK(HUFFMAN_I, 7);
                break;
            case 'a':
                PACK(HUFFMAN_a, 7);
                break;
            case 't':
                PACK(HUFFMAN_t, 7);
                break;
            case 'g':
                PACK(HUFFMAN_g, 7);
                break;
            case 'c':
                PACK(HUFFMAN_c, 7);
                break;
            case 'n':
                PACK(HUFFMAN_n, 7);
                break;
            case 'm':
                PACK(HUFFMAN_m, 7);
                break;
            case 'r':
                PACK(HUFFMAN_r, 7);
                break;
            case 's':
                PACK(HUFFMAN_s, 7);
                break;
            case 'v':
                PACK(HUFFMAN_v, 7);
                break;
            case 'w':
                PACK(HUFFMAN_w, 7);
                break;
            case 'y':
                PACK(HUFFMAN_y, 7);
                break;
            case 'h':
                PACK(HUFFMAN_h, 7);
                break;
            case 'k':
                PACK(HUFFMAN_k, 7);
                break;
            case 'd':
                PACK(HUFFMAN_d, 7);
                break;
            case 'b':
                PACK(HUFFMAN_b, 7);
                break;
            case 'i':
                PACK(HUFFMAN_i, 7);
                break;
            default:
                throw __FILE__ ":deflate_dna_seq: Unknown symbol";
        }
    }
    
    PACK(HUFFMAN_END, 3);

    // Make the bits left-aligned within the terminall byte
    if(bit != 0){

        c = c << (8 - bit);
        ret.push_back(c);
    }

	return ret;
}

// Decompress a DNA sequence
inline bool pop_bit(std::string::const_iterator &m_iter, char &m_offset)
{
    const bool ret = (*m_iter >> m_offset) & 1;
    
    --m_offset;

    if(m_offset < 0){

        m_offset = 7;
        ++m_iter;
    }

    return ret;
}

inline std::string inflate_dna_seq(const std::string &m_bits)
{
	std::string ret;

    if( m_bits.empty() ){
        return std::string();
    }

    std::string::const_iterator i = m_bits.begin();

    // Read bits from left to right
    char offset = 7;

    while(true){

        if(pop_bit(i, offset) == 0){ // 0b0
            if(pop_bit(i, offset) == 0){ // 0b00
                if(pop_bit(i, offset) == 0){ // 0b000
                    ret.push_back('A');
                }
                else{ // 0b001
                    ret.push_back('T');
                }
            }
            else{
                if(pop_bit(i, offset) == 0){ // 0b010
                    ret.push_back('G');
                }
                else{ // 0b011
                    ret.push_back('C');
                }
            }
        }
        else{ // 0b1
            if(pop_bit(i, offset) == 0){ // 0b10
                if(pop_bit(i, offset) == 0){ // 0b100
                    return ret;
                }
                else{ // 0b101
                    if(pop_bit(i, offset) == 0){ // 0b1010
                        if(pop_bit(i, offset) == 0){ // 0b10100
                            if(pop_bit(i, offset) == 0){ // 0b101000
                                ret.push_back('N');
                            }
                            else{ // 0b101001
                                ret.push_back('3');
                            }    
                        }
                        else{ // 0b10101
                            if(pop_bit(i, offset) == 0){ // 0b101010
                                ret.push_back('5');
                            }
                            else{ // 0b101011
                                ret.push_back(' ');
                            }
                        }
                    }
                    else{ // 0b1011
                        if(pop_bit(i, offset) == 0){ // 0b10110
                            if(pop_bit(i, offset) == 0){ // 0b101100
                                ret.push_back('\n');
                            }
                            else{ // 0b101101
                                ret.push_back('\'');
                            }    
                        }
                        else{ // 0b10111
                            if(pop_bit(i, offset) == 0){ // 0b101110
                                ret.push_back('|');
                            }
                            else{ // 0b101111
                                ret.push_back(':');
                            }
                        }
                    }
                }
            }
            else{ // 0b11
                if(pop_bit(i, offset) == 0){ // 0b110
                    if(pop_bit(i, offset) == 0){ // 0b1100
                        if(pop_bit(i, offset) == 0){ // 0b11000
                            if(pop_bit(i, offset) == 0){ // 0b110000
                                if(pop_bit(i, offset) == 0){ // 0b1100000
                                    ret.push_back('-');
                                }
                                else{ // 0b1100001
                                    ret.push_back('M');
                                }
                            }
                            else{ // 0b110001
                                if(pop_bit(i, offset) == 0){ // 0b1100010
                                    ret.push_back('R');
                                }
                                else{ // 0b1100011
                                    ret.push_back('S');
                                }
                            }
                        }
                        else{ // 0b11001
                            if(pop_bit(i, offset) == 0){ // 0b110010
                                if(pop_bit(i, offset) == 0){ // 0b1100100
                                    ret.push_back('V');
                                }
                                else{ // 0b1100101
                                    ret.push_back('W');
                                }
                            }
                            else{ // 0b110011
                                if(pop_bit(i, offset) == 0){ // 0b1100110
                                    ret.push_back('Y');
                                }
                                else{ // 0b1100111
                                    ret.push_back('H');
                                }
                            }
                        }
                    }
                    else{ // 0b1101
                        if(pop_bit(i, offset) == 0){ // 0b11010
                            if(pop_bit(i, offset) == 0){ // 0b110100
                                if(pop_bit(i, offset) == 0){ // 0b1101000
                                    ret.push_back('K');
                                }
                                else{ // 0b1101001
                                    ret.push_back('D');
                                }
                            }
                            else{ // 0b110101
                                if(pop_bit(i, offset) == 0){ // 0b1101010
                                    ret.push_back('B');
                                }
                                else{ // 0b1101011
                                    ret.push_back('I');
                                }
                            }
                        }
                        else{ // 0b11011
                            if(pop_bit(i, offset) == 0){ // 0b110110
                                if(pop_bit(i, offset) == 0){ // 0b1101100
                                    ret.push_back('a');
                                }
                                else{ // 0b1101101
                                    ret.push_back('t');
                                }
                            }
                            else{ // 0b110111
                                if(pop_bit(i, offset) == 0){ // 0b1101110
                                    ret.push_back('g');
                                }
                                else{ // 0b1101111
                                    ret.push_back('c');
                                }
                            }
                        }
                    }
                }
                else{ // 0b111
                    if(pop_bit(i, offset) == 0){ // 0b1110
                        if(pop_bit(i, offset) == 0){ // 0b11100
                            if(pop_bit(i, offset) == 0){ // 0b111000
                                if(pop_bit(i, offset) == 0){ // 0b1110000
                                    ret.push_back('n');
                                }
                                else{ // 0b1110001
                                    ret.push_back('m');
                                }
                            }
                            else{ // 0b111001
                                if(pop_bit(i, offset) == 0){ // 0b1110010
                                    ret.push_back('r');
                                }
                                else{ // 0b1110011
                                    ret.push_back('s');
                                }
                            }
                        }
                        else{ // 0b11101
                            if(pop_bit(i, offset) == 0){ // 0b111010
                                if(pop_bit(i, offset) == 0){ // 0b1110100
                                    ret.push_back('v');
                                }
                                else{ // 0b1110101
                                    ret.push_back('w');
                                }
                            }
                            else{ // 0b111011
                                if(pop_bit(i, offset) == 0){ // 0b1110110
                                    ret.push_back('y');
                                }
                                else{ // 0b1110111
                                    ret.push_back('h');
                                }
                            }
                        }
                    }
                    else{ // 0b1111
                        if(pop_bit(i, offset) == 0){ // 0b11110
                            if(pop_bit(i, offset) == 0){ // 0b111100
                                if(pop_bit(i, offset) == 0){ // 0b1111000
                                    ret.push_back('k');
                                }
                                else{ // 0b1111001
                                    ret.push_back('d');
                                }
                            }
                            else{ // 0b111101
                                if(pop_bit(i, offset) == 0){ // 0b1111010
                                    ret.push_back('b');
                                }
                                else{ // 0b1111011
                                    ret.push_back('i');
                                }
                            }
                        }
                        else{ // 0b11111
                            throw __FILE__ ":inflate_dna_seq: Unexpected symbol";
                        }
                    }
                }
            }
        }
    }

	return ret;
}

#endif // __COMPRESS_DNA
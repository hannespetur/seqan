#ifndef SEQAN_BAM_IO_CRAM_ALIGNMENT_RECORD_H_
#define SEQAN_BAM_IO_CRAM_ALIGNMENT_RECORD_H_

namespace seqan {

class CramSequenceRecord
{
  public:
    String<char> qName;
    String<Iupac> seq;
};


class CramAlignmentRecord
{
  public:
    __uint32 _qID;
    String<CigarElement<> > cigar;
    String<char> qName;
    String<Iupac> seq;
    String<char> qual;
    CharString tags;
    CharString _buffer;

    static __int32 const INVALID_POS = -1;
    static __int32 const INVALID_REF_ID = -1;
    static __int32 const INVALID_LEN = 0;
    static __uint32 const INVALID_QID = 4294967295u;

    // CramAlignmentRecord() : _qID(MaxValue<unsigned>::VALUE) { clear(*this); }
};

} // namespace seqan

#endif // SEQAN_BAM_IO_CRAM_ALIGNMENT_RECORD_H_

#ifndef SEQAN_HTS_IO_HTS_ALIGNMENT_RECORD_H_
#define SEQAN_HTS_IO_HTS_ALIGNMENT_RECORD_H_

namespace seqan {

class HtsSequenceRecord
{
  public:
    String<char> qName;
    String<Iupac> seq;
};


class HtsAlignmentRecord
{
  public:
    __uint32 _qID;
    String<CigarElement<> > cigar;
    String<char> qName;
    String<Iupac> seq;
    String<char> qual;
    String<char> tags;
    String<char> _buffer;

    static __int32 const INVALID_POS = -1;
    static __int32 const INVALID_REF_ID = -1;
    static __int32 const INVALID_LEN = 0;
    static __uint32 const INVALID_QID = 4294967295u;

    // CramAlignmentRecord() : _qID(MaxValue<unsigned>::VALUE) { clear(*this); }
};

} // namespace seqan

#endif // SEQAN_HTS_IO_HTS_ALIGNMENT_RECORD_H_

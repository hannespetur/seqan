#ifndef SEQAN_HTS_IO_HTS_FILE_H_
#define SEQAN_HTS_IO_HTS_FILE_H_

#include <seqan/hts_io/hts_alignment_record.h>

#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>


namespace seqan {

class HtsFile
{
  public:
    const char * filename;
    htsFile * fp;
    bam_hdr_t * hdr;
    bam1_t * hts_record;
    hts_idx_t * hts_index;
    hts_itr_t * hts_iter = NULL;

    HtsFile(const char * filename)
        : filename(filename)
        {}
};

inline void
_parseHtsSequenceRecord(HtsSequenceRecord & record, HtsFile & file)
{
    record.qName = bam_get_qname(file.hts_record);
    int32_t lqseq = file.hts_record->core.l_qseq;
    resize(record.seq, lqseq);
    uint8_t* seqptr = bam_get_seq(file.hts_record);

    for (int i = 0; i < lqseq; ++i)
    {
      record.seq[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
    }
}

inline bool
readRecord(HtsSequenceRecord & record, HtsFile & file)
{
    if (sam_read1(file.fp, file.hdr, file.hts_record) < 0)
        return false;

    _parseHtsSequenceRecord(record, file);
    return true;
}

inline void
open(HtsFile & file)
{
    file.fp = hts_open(file.filename, "r");
    if (file.fp == NULL)
    {
        SEQAN_FAIL("Could not open file with filename %s", file.filename);
    }
    
    file.hdr = sam_hdr_read(file.fp);
    file.hts_record = bam_init1();
}

inline bool
loadIndex(HtsFile & file)
{
    file.hts_index = sam_index_load(file.fp, file.filename);
    return file.hts_index != NULL;
}

inline bool
loadIndex(HtsFile & file, bool build_if_unavailable)
{
    if (build_if_unavailable)
    {
        file.hts_index = sam_index_load(file.fp, file.filename);

        if (file.hts_index == NULL)
        {
            bam_index_build(file.filename, 1000);
            return loadIndex(file);
        }

        return true;
    }

    return loadIndex(file);
}

inline void
setRegion(HtsFile & file, const char * region)
{
    if (file.hts_iter != NULL)
        sam_itr_destroy(file.hts_iter);

    file.hts_iter = sam_itr_querys(file.hts_index, file.hdr, region);
}

inline bool
readRegion(HtsSequenceRecord & record, HtsFile & file)
{
    if (sam_itr_next(file.fp, file.hts_iter, file.hts_record) < 0)
        return false;

    _parseHtsSequenceRecord(record, file);
    return true;
}

} // namespace seqan


#endif  // SEQAN_HTS_IO_HTS_FILE_H_

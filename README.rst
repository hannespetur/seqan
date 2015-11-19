SeqAn - The Library for Sequence Analysis (now with HTS!)
=========================================================

This is a fork of SeqAn, with added extensibility of the htslib. List of new added features:
 * BCF file support.
    - Reading implemented.
    - Writing **not** implemented.
    - Indexes (Tabix)
        + Reading implemented.
        + Building **not** implemented.
    - Parsing to SeqAn's VCF record is implemented.
 * HTS (SAM/BAM/CRAM) support.
    - Reading implemented.
    - Writing implemented.
    - Indexes (.bai, .csi, .crai)
        + Reading implemented.
        + Building implemented.
    - Parsing to SeqAn alignment records.
        + Parsing of qName, seq, qual and cigar implemented. The rest is **not**.


To make use of these features you need to add `-DSEQAN_USE_HTSLIB=1` to your CXX compiler flags. If you need a feature which is not implemented you can ask me (Hannes Pétur, hannese@decode.is) to add it, or use the `htslib API <https://github.com/samtools/htslib>`_ directly.


Examples
--------------
BCF/Tabix example
~~~~~~~~~~~~~~~~~
.. code-block:: cpp

  seqan::Tabix index;
  seqan::open(index, "/path/to/my/file/example.vcf.gz");
  seqan::setRegion(index, "chrX:A-B"); // "chrX" and "chrX:A" also supported
  seqan::VcfRecord record;
  
  while (seqan::readRegion(record, index))
  {
    // Do stuff with record
  }


HTS file read and write example (SAM/BAM/CRAM is automatically detected)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: cpp

  seqan::HtsFile read_file("/path/to/some/existing/file.cram", "r");
  seqan::HtsFile write_file("/path/to/a/new/file.cram", "wb"); // binary mode required for BAM and CRAM

  seqan::copyHeader(write_file, read_file);
  seqan::writeHeader(write_file);

  while (seqan::readRecord(read_file))
  {
    seqan::copyRecord(write_file, read_file);
    // Here you could change the write_file.hts_record if you want.
    seqan::writeRecord(write_file);
  }
  // Here we have copied all records of read_file and written them to write_file


HTS file index example
~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: cpp

  seqan::HtsFile hts_file("/path/to/some/existing/file.cram", "r");

  if (!seqan::loadIndex(hts_file))
  {
    // Build it if we cannot find it
    seqan::buildIndex(hts_file);
    seqan::loadIndex(hts_file);
  }

  seqan::setRegion(hts_file, "chrX:A-B");
  seqan::HtsSequenceRecord record; // Only parses qName and sequence, use seqan::HtsAlignmentRecord to parse all

  while (readRegion(record, hts_file))
  {
    // Do stuff with each record that overlaps the chrX:A-B region.
  }


SeqAn
--------------

For information about SeqAn check out `<https://github.com/seqan/seqan>`_


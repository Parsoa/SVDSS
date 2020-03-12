# Variation-Discovery

Some (maybe evident) notes:
* we must use some sort of index
* before counting the number of solutions, I would like to see how many reads are effectively different between parents and child
* we can try to reduce the number of reads and then analyze the reads we are sure are "different". I expect to have some (maybe a lot?) reads that are shared between the parents and the child (but I'm not sure - I haven't worked with HiFi data and SVs so I don't know what to expect)
* a preliminary (not perfect) analysis can help us in understand what we have to do

**First approach:**
1. we index the parents' read with a suffix tree (as suggested by Rayan)
2. we search for the reads from the child and we keep only those reads that are not completely in the tree

**Second approach (I'm not sure this is feasible):**
1. we create the references of the parents by combining the reference genome and the VCF containing the known SVs (with SNPs and indels this can be done with bcftools. We have to check if this works also for SVs)
2. we align the reads from the child to the references. I expect to have clipped reads near a SV (maybe we can just use the provided BAMs?)
3. we analyze the clipped reads to extract unique reads

This second approach may be too biased towards the known set of SVs but I think it could give us an rough idea of what to expect.

---

Both these approaches:
* assume that some reads are shared among the parents and the child
* require a post processing step to detect unique strings (we have to manage sequencing errors somehow)

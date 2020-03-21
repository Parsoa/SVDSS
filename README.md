# Variation-Discovery

Some (maybe evident) notes:
* we must use some sort of index
* before counting the number of solutions, I would like to see how many reads are effectively different between parents and child
* we can try to reduce the number of reads and then analyze the reads we are sure are "different". I expect to have some (maybe a lot?) reads that are shared between the parents and the child (but I'm not sure - I haven't worked with HiFi data and SVs so I don't know what to expect)
* a preliminary (not perfect) analysis can help us in understand what we have to do

**First approach:**
1. we index the parents' read with a suffix tree (as suggested by Rayan)
2. we search for the reads from the child and we keep only those reads that are not completely in the tree
    * Question here is what to insert to the tree from the child reads. Assume we have read P in the parent and read R in the child and R[1:i] is found in P. We can just try matching R[i + j:] with P until we find a j that matches with length of at least 200bp and then insert R[i:j] into the tree. Problem is R[i:j] might just be too short to be meaningful and some sequence around it also needs to be considered.
        * How to find how much sequence surrounding R[i:j] to consider?
            * possible solution: change the granularity of match, say always match in chunks of 8bp each. A mismatch will be at least 16bp long which is more likely to already not exist in parent. (generally 32bp should be a better lowerbound)

**Approach 1.1:**
1. we index the parents' samples with an FM-index
2. we perform kmer counting on the child sample (for different values of k, let's say from 17 to 200/250? For small values of k we can use KMC. I have to check if we can use KMC for bigger values of k)
3. we search each kmer in the index and we count how many kmers are only in the child
4. we can plot a sort-of histogram to see the distribution of unique kmers

My hope is that we can use such unique kmers to "assembly" longer unique strings.

**Second approach (I'm not sure this is feasible - EDIT 03/21: I don't think it's feasible. I'll try approach 1.1 first):**
1. we create the references of the parents by combining the reference genome and the VCF containing the known SVs (with SNPs and indels this can be done with bcftools. We have to check if this works also for SVs)
2. we align the reads from the child to the references. I expect to have clipped reads near a SV (maybe we can just use the provided BAMs?)
    * If there is a SV in the child (compared to ref), there will be clipped reads in the BAM that map near the breakpoints. Maybe use those?
3. we analyze the clipped reads to extract unique reads

This second approach may be too biased towards the known set of SVs but I think it could give us an rough idea of what to expect.

---

Both these approaches:
* assume that some reads are shared among the parents and the child
* require a post processing step to detect unique strings (we have to manage sequencing errors somehow)

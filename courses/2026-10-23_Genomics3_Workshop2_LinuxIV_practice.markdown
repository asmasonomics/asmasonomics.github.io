---
layout: page
permalink: /courses/Genomics3_Workshop2_LinuxIV_practice_Oct2026
---

![Genomics3 banner](/assets/coursefiles/2023_Genomics/genomics_banner.jpeg){:class="img-responsive"}

<span style="font-size:1.6em;">**Genomics 3 - Workshop 2: LinuxIV -- Practice Exercise**</span><br/>

<p align="justify">
Well done for completing the taught material from Workshop 2. We've written some exercises here so you can revise and check your understanding. At the end, you can run a little script which will check your answers.<br/><br/>
Ask us if you're not sure. Take some time to explore and get used to working in Linux - you will use it in each data workshop and for your assessment. Students regularly come unstuck with typos and by getting confused about which diretory they are in (or their data is in). Practice will increase your confidence!<br/>
</p>

## Practice material
<p align="justify">
Work in your own terminal and keep a lab book record as you go.<br/><br/>
You will be using skills from the earlier workshop material: navigating directories, <code>pwd</code>, symbolic links, wildcards, <code>grep</code>, redirection into files with <code>></code> and pipes <code>|</code>.<br/><br/>
You'll also meet one new command, <code>zcat</code>, and learn how to do basic maths at the command line.
</p>

### Setup 

1\. Make sure you're in your own directory

```sh
cd /shared/biology/bioldata1/bl-00087h/students/$USER
pwd
```

2\. Make a new directory called <code>practice_w2</code> and move into it.

### Part A - Link the data

3\. Create symbolic "soft" links to all the files in the *L. braziliensis* genome directory (`/shared/biology/bioldata1/bl-00087h/data/L.braziliensis/genome/`) where the filename starts with Lbraz.
You did something very similar to this earlier.

4\. Run `ls -lF` to confirm you now have a `.fasta` file and at least one `.fastq.gz` file linked into `practice_w2`, each with an arrow (`->`) pointing back to the original.

5\. We only care about the files which end in `.fasta` or `.fastq.gz`. Remove the rest in one command. Can you work out how you could have created the soft links only to these files in step 3?

### Part B - Count the FASTA sequences
6\. Count how many sequences are in the `.fasta` file. Remember, in FASTA format each sequence name (header) is on a line starting with `>`.

7\. Save that number to a new file called `fasta_seq_count.txt` (don't just print it to the screen — redirect it into the file).

### Part C - Look at a FASTQ file, then subset
FASTQ files are usually gzip-compressed (`.gz`), so you can't just `less` or `head` them directly and expect to read them. There's a command called `zcat` that decompresses a `.gz` file "on the fly" and streams its contents straight to the screen (or into a pipe) — without ever writing an unzipped copy to disk.

8\. Look at the first few lines of one of your linked `.fastq.gz` files, using `zcat` piped into a command you already know:

```sh
zcat Lbraz-subset.fastq.gz | head
```

Look carefully at the structure. Every read in a FASTQ file is written as **4 lines**: a header line (starting with `@`), the sequence, a `+` line (sometimes a repeat of the read name), and a quality score line.<br/>
This is easier to see using the `less` command followed by the `-S` flag. This allows lines to run off the screen rather than wrapping onto the next line. Press `q` to escape.

```sh
zcat Lbraz-subset.fastq.gz | less -S
```

9\. Now, work out how many **reads** are in this fastq file (not how many lines!). Save that number to a file called `read_count.txt`. Hint - how does number of lines relate to number of reads in a FASTQ file?

10\. Check your answer! Is this number correct? Assuming you use `Lbraz-subset.fastq.gz`, do the following:

```sh
zcat Lbraz-subset.fastq.gz | wc -l
```

This should give an answer of `113740`. <br/>
Now, let's do some basic command line maths. It looks a bit weird, but it follows `$(( ))` nomenclature, where the maths goes in the middle of the double brackets. This is called **arithmetic expansion** in linux.

```sh
echo $(( 113740 / 4 ))
```

Is this the same answer as in `read_count.txt`? No.<br/><br/>
As we'll cover in later workshops, the "quality" of a base (i.e. confidence the sequencer got it right) is encoded as an ASCII character, including the `@`. By chance, the first base may have the quality `@` meaning you count lines which are not just headers. The `grep` approach only works if you can be certain your pattern matching is unique.<br/><br/>
But now we have the maths, and we can also use <code>` `</code>. In coding the different quotes have different meanings - <code>' " `</code> are different. In linux <code>` `</code> is used to say "do everything inside these back ticks first", so we can do the following:

```sh
echo $(( `zcat Lbraz-subset.fastq.gz | wc -l` / 4 )) > read_count.txt
```

This should now match up! Even when we are confident with a command, or we've done something before it is **so so so** important to double check. We can't necessarily always imagine every scenario, every edge case, every weird thing that someone else has done to our dataset. Always check. It will save you lots of time in the future.

11\. Create a new file called `first10reads.fastq` that contains **only** the first 10 reads of the fastq file.
Also see that now the file is `fastq` not `fastq.gz`, so we're back into plain text so you can use cat rather than `zcat`.

12\. Check your work: count how many lines in `first10reads.fastq` start with the `@` symbol (the way you'd check for `>` in a FASTA file), and save that count to a file called `header_count.txt`.<br/>
Is the number what you expected? If not — have a look at the contents of `first10reads.fastq` with `less -S` and think about the maths we did above.

### Check your answers
Once you've done all of the above, run the checking script from inside your `practice_w2` directory:

```sh
sh /shared/biology/bioldata1/bl-00087h/data/check_workshop2_practice.sh
```

It will tell you `[PASS]` or `[FAIL]` for each part, and give you a hint if something's wrong. You don't need to get everything right first time. The key thing is that you understand if you have gone wrong and know how to fix it.

## Finishing up
<p align="justify">
Assuming you got a set of PASS values from the script - well done! You've completed the Linux training and are now ready for the data workshops over the next few weeks.<br/><br/>
Remember, the DataCamp access you have is unlimited until the end of January, so you can continue to develop your coding skills with certificate evidence for your CV if you want to. After January your account will remain but your access will become more limited. 
</p>

---
layout: page
permalink: /courses/MasonLab_Linux_onboarding
---

<p align="justify">
Welcome to the group!<br/><br/>
Much of your work in the lab will involve using the University's high performance computing infrastructure: Viking2. The is a Linux system. You may be familiar with Linux already. If so, this page should just be a series of signposts for you. If not, you get signposts to training as well.
</p>

### Useful links
<a href="https://asmasonomics.github.io/">Group website - asmasonomics.github.io</a><br/>
<a href="https://vikingdocs.york.ac.uk/index.html">Viking2 documentation</a><br/>

### Viking access form and existing code
<p align="justify">
To access Viking2, please complete <a href="https://docs.google.com/forms/d/e/1FAIpQLSfXkL10ypU6EQCBB2jS5oDwTpRMo77ppl7dvdbLnXm5zrKR7Q/viewform">the Viking2 application form</a>. For this you will need a <code>project code</code> - <a class="u-email" href="mailto:andrew.mason@york.ac.uk">ask Andrew for this.</a><br/><br/>
As default, each user is allocated 2TB of storage space in the <code>/mnt/scratch/users/your-UUN/</code> directory. As soon as your user access is approved, <a class="u-email" href="mailto:itservices@york.ac.uk">please email IT services</a> to increase your quota. Use the subject heading "Viking - quota increase for [your-UUN]" and ask to increase your quota to 10TB unless Andrew has said you should ask for more.<br/><br/>
You should also <a href="https://github.com/">set up a GitHub account</a> if you don't already have one. As academics often move around a lot, it is prudent to link this to a personal email address as well as your York one, and to make your username something you like. Once you have this, let Andrew know and he will add you to our <a href="https://github.com/Mason-Lab-Code">Mason Lab GitHub organisation</a>, where you can access code and link any project repositories. On here, Richard has put some very useful information about how we organise our projects, and also some links to software carpentry introductory courses - use these!<br/><br/>
If you're using a managed Windows machine, and you're here longer than a Summer project, <a class="u-email" href="mailto:david.nelmes@york.ac.uk">please email David Nelmes</a>, cc'ing in Andrew, and ask if you can have elevated admin rights on your machine as you're a bioinformatician needing to install lots of new software outside the SoftwareCentre.
<br/>
</p>

### Accessing Viking
<p align="justify">
If you're on campus, great. If not, you will need to use the <a href="https://www.york.ac.uk/it-services/tools/vpn/">VPN service</a> with 2FA.<br/><br/>
From Windows you can use the built in PowerShell or  <a href="https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html">download the latest 64-bit puTTY SSH client.</a><br/>
If you're on a Mac use the Terminal, and if on a Linux machine you won't be needing this tutorial.<br/><br/>
It is prudent to also have some kind of SCP/SFTP GUI client to move files between Viking and your desktop. For Windows the best one is <a href="https://winscp.net/eng/download.php">WinSCP</a>, but <a href="https://filezilla-project.org/download.php?type=client">FileZilla Client</a> (not FileZilla Server...) is OK too (and the best one for Mac). <br/><br/>
When you log on to Viking, you will be in <code>/users/your-UUN/</code>. Don't do any work here. It has almost no storage. You should immediately change directory to either your own personal scratch space <code>/mnt/scratch/users/your-UUN/</code> or to the group space <code>/mnt/scratch/projects/biol-cancerinf-2020/</code><br/>
To save your fingers, we recommend you create symbolic links:
<br/>
</p>

```sh 

# long way to change into the group space
cd /mnt/scratch/projects/biol-cancerinf-2020/

# instead, create symbolic links
cd 
ln -s /mnt/scratch/projects/biol-cancerinf-2020/ group
ln -s /mnt/scratch/users/your-UUN/ scratch

```

<p align="justify">
The only files you will likely ever edit in your initial login space <code>/users/your-UUN/</code> will be <code>~/.bash_profile</code> where you can add <code>module load</code> statements, and <code>~/.bashrc</code> where you may add environment variables, including aliases. <br/><br/>
If none of that is clear to you, don't worry (and don't mess with them until you do).
<br/>
</p>

### Key messages
1. Refer to the documentation
2. Have good data management - refer to Richard's style guides
3. Work within scratch (personal or group)
4. Remember, there is no recycle bin - if you delete stuff <b>it is gone forever</b>.

### Next steps
<p align="justify">
If you have worked on Linux before and are feeling good, read the Viking documentation (selectively) and get cracking! Pay particular attention to available modules and how the <code>slurm</code> scheduler works if you haven't used that before.<br/><br/>
Always ask Andrew and Richard for help if you're unsure - this is not a problem. Far better to ask early and get the help you need rather than struggle for ages or worsen issues you might have.<br/><br/>
If you've not worked with Linux before, don't worry - we have some training and guidance for you. As with any coding language it can take a bit of getting used to. But you'll get there.
<br/>
</p>

### Linux terminal keyboard shortcuts
<code>Ctrl A</code> - return to start of terminal line<br/>
<code>Ctrl E</code> - go to end of terminal line<br/>
<code>Ctrl U</code> - delete from cursor location to start of code line<br/>
<code>Ctrl R</code> - quick search tool through last code entries<br/>
<code>Tab</code> - autocomplete, double tap to see all possible options<br/>
<kbd>&uarr;</kbd> - scroll back up through previous commands<br/>

### Coding training
<p align="justify">
As mentioned, Richard has added some links on the <a href="https://github.com/Mason-Lab-Code">Mason Lab GitHub organisation</a>. We have also set up a DataCamp subscription for the group - ask Andrew for access - <i>2024-04-05 awaiting confirmation - more details to follow...</i><br/><br/>
For broader training courses, check out <a href="https://training.galaxyproject.org/">the Galaxy Training Network (GTN)</a> and <a href="https://tess.elixir-europe.org/">Elixir's TeSS</a>.<br/>
</p>
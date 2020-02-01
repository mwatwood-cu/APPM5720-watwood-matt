# APPM 5720 Matt Watwood Project
## Intro Fun
Haiku's are easy<br>
But they don't always make sense<br>
Refrigerator<br><br>

## Description
I will be working on implementing the embedded boundary layer method listed for those of us without a direct research project. <br><br>

It seems to me that given the effect of nearest neighbors on the boundary method, it can be parallelized by separating the domain and using ghost points along the edges including the boundary. This will make the entire domain separable and able to be processed on separate processors, then only a small amount of node to node communication is needed to transfer the ghost points. I have not dug into the provided paper enough to guarentee this method will work, but as a first glance and read through I think this will work.<br><br>

Realistically, I would plan to spend 1-4 hours per week on this project on an average week. Starting familiarizing myself with the paper, then implementing it in C++, and finally looking to parallelize.



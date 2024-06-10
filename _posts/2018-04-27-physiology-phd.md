---
layout: post
title: "Beats of Stress"
date: 2018-04-27 21:01:00
description: "A physiological perspective on my PhD exams"
tags: PhD-journey, data-visualization
categories: miscellany
thumbnail: assets/img/posts/comps_hr.png
giscus_comments: true
---

The [Comprehensive Exam](https://en.wikipedia.org/wiki/Comprehensive_examination) (more affectionately known as "comps") is a notoriously stressful milestone in the PhD journey. It is a high-stakes exam that tests the student's knowledge of their field and their ability to think critically. The format of the exam varies by departments, but usually consists of both a written and an oral component. For my doctoral program, this involved writing a research proposal laying out the plans for my dissertation research, and then defending it in front of my doctoral committee, which comprised of five faculty members from related fields.

As a competitive cyclist during my PhD, I regularly used heart rate metrics to track my training and recovery. As exam anxiety started mounting in the days leading up to the exam, my curiosity piqued about how my body would respond physiologically to this gruelling ordeal. To explore this, I decided to wear my heart rate monitor to the exam to gather some data.

Now, onto the data: the physiological profile of my heart was quite revealing. As the exam commenced, my heart rate spiked sharply to 140 bpm, well within the fat-burning zone as I stood before my doctoral committee with my presentation behind me. This was a clear sign of the pre-exam jitters. The anticipation of the questions, the intimidating presence of the professors, the sheer weight of the moment - all acting as conductors in this orchestra of anxiety.

<div class="row justify-content-center mt-3">
    <div class="col-12 mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/posts/comps_hr.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
<div class="caption">
    Beat by beat, the heart tells the story of my PhD exam.
</div>
</div>

Fortunately, my professors proved to be a compassionate and supportive group, helping to ease my nerves. As I answered their challenging questions with accuracy (or at least confidence), I felt my heart begin to ease, reflected in the gradual decline of my heart rate. Yet, there were still notable spikes, likely triggered by especially tough inquiries, sending my pulse racing once more.

As my heart rate stabilized towards the end of the exam, the ultimate trial loomed: the verdict. They ask you to step outside the room, while the examiners debate your fate. My heart rate surged dramatically again, doing a pole vault during those agonizing moments spent in anticipation outside the room.

In the end, I survived the comps, my dissertation proposal was met with positive enthusiasm and the examiners' insights helped sharpen and refine the questions I was going to explore in my doctoral work. In retrospect, the PhD oral exam, viewed through the lens of a heart rate monitor revealed the physiological response of the body to this intellectual pressure. This curious experiment was a welcome distraction, and though I wouldn't advocate it as a fat-burning regime, it undoubtedly added a unique dimension to my academic odyssey.

<ins>Notes:</ins>

- My resting heart rate during this period was around 45 bpm, while my average heart rate during the exam was 84 bpm, a significant increase.
- Used the [Wahoo TICKR](https://uk.wahoofitness.com/devices/running/heart-rate-monitors/tickr-buy) heart rate monitor to capture the data, exported using the [Wahoo app](https://play.google.com/store/apps/details?id=com.wahoofitness.fitness&hl=en_US), cleaned up the data and plotted the graph in [R](https://www.r-project.org/). Code for the plot available on my [GitHub](https://github.com/abhilesh/Miscellaneous_scripts/tree/master/PhD_Comps_HR).

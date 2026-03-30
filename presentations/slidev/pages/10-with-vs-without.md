---
layout: center
---

# WITH vs WITHOUT Mechanistic Offset

<div class="flex justify-center">
<img src="/images/anim_with_vs_without.gif" class="h-72 rounded shadow-lg" />
</div>

<div class="mt-2 text-xs text-center max-w-2xl mx-auto">

**What you're seeing:** The blue line (WITH I\* offset) tracks seasonal dynamics because the mechanistic prediction provides the transmission shape. The red line (I\*=0, standard geostatistical) collapses to a flat mean — without mechanistic information, it cannot capture temporal structure. This is the core comparison: with mechanistic information vs without (setting I\* = 0).

</div>

<!--
This is the key comparison. The WITH offset model uses I-star as a log-offset, so the GP only needs to capture departures from the mechanistic prediction. The WITHOUT model (I-star equals zero) is the standard geostatistical approach — equivalent to what MAP does without vector data. The GP must learn all spatial and temporal structure from data alone, requiring much more GP variance and producing worse parameter estimates.
-->

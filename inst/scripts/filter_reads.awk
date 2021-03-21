#!/usr/bin/awk -f
($9 != isize && $5 != q) || ($9 == isize && $7 != "=" && $5 != q)

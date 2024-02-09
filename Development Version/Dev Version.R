### Development version of Identify Halogenated Spectra
# source('Dev Version.R', echo = FALSE)

# Put the ChromaTOF MSP files (renamed with .csv extensions) in a folder named 'Input Spectra'. The script will run on every file in that folder.

rm(list = ls())

# Theoretical isotopic distributions for Br 1 to 10, Cl 1 to 10, and combinations of Br and Cl up to Br5Cl5. Colums are theoretical_mz, theoretical_intensity, theoretical_percent_intensity, distribution_name, and peak_number.
# To generate from theoretical_distributions.csv:
# x <- read.csv(file = 'theoretical distributions.csv', stringsAsFactors = FALSE)
# to_remove <- c('theoretical_intensity', 'peak_number')
# x <- x[, !(names(x) %in% to_remove), drop = F]
# dput(x)
t_distributions <- # theoretical distributions
structure(list(theoretical_mz = c(NA, NA, NA, NA, NA, NA, NA, 
NA, NA, NA, 234, 235.01, 236, 237.01, 238.01, NA, NA, NA, NA,
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 311.91, 312.92, 
313.91, 314.92, 315.91, 316.91, 317.92, NA, NA, NA, NA, NA, NA,
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 389.83, 390.83, 391.82,
392.83, 393.82, 394.82, 395.82, 396.82, 397.83, 398.83, NA, NA,
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 467.74, 468.74, 469.73,
470.74, 471.73, 472.74, 473.73, 474.73, 475.73, 476.73, 477.73,
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 545.65, 546.65,
547.64, 548.65, 549.64, 550.65, 551.64, 552.64, 553.64, 554.64,
555.64, 556.64, 557.64, NA, NA, NA, NA, NA, NA, NA, NA, NA, 623.56,
624.56, 625.55, 626.56, 627.55, 628.56, 629.55, 630.55, 631.55,
632.55, 633.55, 634.55, 635.55, 636.55, 637.55, NA, NA, NA, NA, 
NA, NA, NA, 701.47, 702.47, 703.47, 704.47, 705.46, 706.47, 707.46,
708.46, 709.46, 710.46, 711.46, 712.46, 713.46, 714.46, 715.45,
716.46, NA, NA, NA, NA, NA, NA, NA, NA, 779.38, 780.38, 781.38,
782.38, 783.37, 784.38, 785.37, 786.38, 787.37, 788.37, 789.37,
790.37, 791.37, 792.37, 793.36, 794.37, 795.36, 796.36, NA, NA,
NA, NA, NA, NA, 857.29, 858.29, 859.29, 860.29, 861.28, 862.29, 
863.28, 864.29, 865.28, 866.28, 867.28, 868.28, 869.28, 870.28,
871.27, 872.28, 873.27, 874.28, 875.27, 877.28, NA, NA, 935.2,
936.2, 937.2, 938.2, 939.19, 940.2, 941.19, 942.2, 943.19, 944.19, 
945.19, 946.19, 947.19, 948.19, 949.18, 950.19, 951.18, 952.19,
953.18, 954.18, 955.18, 956.18, NA, NA, NA, NA, NA, NA, NA, NA,
NA, NA, NA, NA, 190.05, 191.06, 192.05, 193.06, 194.06, NA, NA,
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
NA, 224.02, 225.02, 226.01, 227.02, 228.01, 229.01, 230.02, 231.02,
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
257.98, 258.98, 259.97, 260.98, 261.97, 262.97, 263.97, 264.97,
265.97, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 291.94,
292.94, 293.94, 294.94, 295.93, 296.94, 297.93, 298.93, 299.93,
300.93, 301.93, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
NA, 325.9, 326.9, 327.9, 328.9, 329.89, 330.9, 331.89, 332.89,
333.89, 334.89, 335.88, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
NA, NA, NA, 359.86, 360.86, 361.86, 362.86, 363.85, 364.86, 365.85,
366.85, 367.85, 368.85, 369.85, 370.85, 371.84, NA, NA, NA, NA, 
NA, NA, NA, NA, NA, NA, NA, 393.82, 394.82, 395.82, 396.82, 397.82,
398.82, 399.81, 400.82, 401.81, 402.81, 403.81, 404.81, 405.8,
406.82, NA, NA, NA, NA, NA, NA, NA, NA, 427.78, 428.79, 429.78,
430.78, 431.78, 432.78, 433.77, 434.78, 435.77, 436.77, 437.77,
438.77, 439.77, 440.77, 441.76, NA, NA, NA, NA, NA, NA, NA, NA,
NA, 461.74, 462.75, 463.74, 464.74, 465.74, 466.74, 467.73, 468.74,
469.73, 470.73, 471.73, 472.73, 473.73, 474.73, 475.72, NA, NA,
NA, NA, NA, NA, NA, NA, NA, 495.7, 496.71, 497.7, 498.7, 499.7,
500.7, 501.7, 502.7, 503.69, 504.7, 505.69, 506.69, 507.69, 508.69,
509.68, 510.69, 511.68, 513.68, NA, NA, NA, NA, NA, NA, NA, NA,
267.97, 268.97, 269.96, 270.97, 271.96, 272.96, 273.97, 274.97,
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
301.93, 302.93, 303.92, 304.93, 305.92, 306.92, 307.92, 308.92, 
309.93, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
NA, 335.89, 336.89, 337.89, 338.89, 339.88, 340.89, 341.88, 342.88,
343.88, 344.88, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
NA, NA, 369.85, 370.85, 371.85, 372.85, 373.84, 374.85, 375.84,
376.84, 377.84, 378.84, 379.84, 380.84, NA, NA, NA, NA, NA, NA,
NA, NA, NA, NA, 403.81, 404.81, 405.81, 406.81, 407.8, 408.81,
409.8, 410.8, 411.8, 412.8, 413.8, 414.8, 415.79, 416.8, NA,
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 345.88, 346.88, 347.87,
348.88, 349.87, 350.87, 351.87, 352.87, 353.88, NA, NA, NA, NA,
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 379.84, 380.84, 381.83,
382.84, 383.83, 384.84, 385.83, 386.83, 387.83, 388.83, 389.83,
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 413.8, 414.8, 415.8, 
416.8, 417.79, 418.8, 419.79, 420.79, 421.79, 422.79, 423.79,
424.79, 425.79, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 447.76,
448.76, 449.76, 450.76, 451.75, 452.76, 453.75, 454.75, 455.75,
456.75, 457.75, 458.75, 459.74, 461.75, NA, NA, NA, NA, NA, NA,
NA, NA, NA, NA, 481.72, 482.72, 483.72, 484.72, 485.72, 486.72,
487.71, 488.72, 489.71, 490.71, 491.71, 492.71, 493.7, 494.71,
495.7, NA, NA, NA, NA, NA, NA, NA, NA, NA, 423.79, 424.79, 425.78,
426.79, 427.78, 428.79, 429.78, 430.78, 431.78, 432.78, 433.78,
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 457.75, 458.75,
459.75, 460.75, 461.74, 462.75, 463.74, 464.74, 465.74, 466.74,
467.74, 468.74, 469.74, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
NA, 491.71, 492.71, 493.71, 494.71, 495.7, 496.71, 497.7, 498.7,
499.7, 500.7, 501.7, 502.7, 503.69, 504.7, NA, NA, NA, NA, NA,
NA, NA, NA, NA, NA, 525.67, 526.67, 527.67, 528.67, 529.66, 530.67,
531.66, 532.67, 533.66, 534.66, 535.66, 536.66, 537.65, 538.66,
539.65, NA, NA, NA, NA, NA, NA, NA, 559.63, 560.63, 561.63, 562.63,
563.63, 564.63, 565.62, 566.63, 567.62, 568.62, 569.62, 570.62, 
571.62, 572.62, 573.61, 574.63, 575.61, NA, NA, NA, NA, NA, NA,
NA, NA, NA, 501.7, 502.7, 503.69, 504.7, 505.69, 506.7, 507.69,
508.69, 509.69, 510.69, 511.69, 512.69, 513.69, NA, NA, NA, NA,
NA, NA, NA, NA, NA, NA, NA, 535.66, 536.66, 537.66, 538.66, 539.65,
540.66, 541.65, 542.65, 543.65, 544.65, 545.65, 546.65, 547.64,
548.65, NA, NA, NA, NA, NA, NA, NA, NA, 569.62, 570.62, 571.62,
572.62, 573.61, 574.62, 575.61, 576.62, 577.61, 578.61, 579.61,
580.61, 581.61, 582.61, 583.6, 584.61, NA, NA, NA, NA, NA, NA,
NA, NA, 603.58, 604.58, 605.58, 606.58, 607.58, 608.58, 609.57,
610.58, 611.57, 612.57, 613.57, 614.57, 615.57, 616.57, 617.56,
618.57, 619.56, 620.56, NA, NA, NA, NA, NA, NA, 637.54, 638.54,
639.54, 640.54, 641.54, 642.54, 643.53, 644.54, 645.53, 646.53,
647.53, 648.53, 649.53, 650.53, 651.52, 652.53, 653.52, 654.52,
NA, NA, NA, NA, NA, NA, 579.61, 580.61, 581.61, 582.61, 583.6,
584.61, 585.6, 586.6, 587.6, 588.6, 589.6, 590.6, 591.6, 592.6,
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 613.57, 614.57, 615.57, 
616.57, 617.56, 618.57, 619.56, 620.57, 621.56, 622.56, 623.56,
624.56, 625.56, 626.56, 627.55, 628.56, 629.56, NA, NA, NA, NA,
NA, NA, NA, 647.53, 648.53, 649.53, 650.53, 651.53, 652.53, 653.52,
654.53, 655.52, 656.52, 657.52, 658.52, 659.52, 660.52, 661.51,
662.52, 663.51, NA, NA, NA, NA, NA, NA, NA, 681.49, 682.49, 683.49,
684.49, 685.49, 686.49, 687.48, 688.49, 689.48, 690.48, 691.48,
692.48, 693.48, 694.48, 695.47, 696.48, 697.47, 698.47, NA, NA,
NA, NA, 715.45, 716.45, 717.45, 718.45, 719.45, 720.45, 721.44,
722.45, 723.44, 724.45, 725.44, 726.44, 727.44, 728.44, 729.44,
730.44, 731.43, 732.43, 733.43, 734.43, NA, NA), theoretical_percent_intensity = c(0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 13.37, 96.92, 12.34, 0.58, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 49.95, 7.23,
100, 13.24, 49.8, 6.12, 0.41, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 34.68, 4.92, 100, 12.76, 98.83, 12.76, 32.04,
3.99, 0.27, 0.03, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 18.59,
2.3, 67.12, 8.82, 100, 13.06, 67.06, 9.25, 14.05, 1.87, 0.09, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10.24, 1.14, 50.55, 5.94,
100, 12.69, 97.97, 11.62, 48.1, 6.65, 9.24, 1.21, 0.14, 0, 0,
0, 0, 0, 0, 0, 0, 0, 4.46, 0.83, 31.58, 4.38, 77.25, 9.59, 100,
12.9, 71.51, 9.59, 27.81, 3.92, 4.99, 0.47, 0.04, 0, 0, 0, 0,
0, 0, 0, 2.68, 0.28, 19.92, 2.72, 62.56, 6.94, 98.22, 11.56,
100, 10.59, 58.3, 8.4, 18.38, 1.78, 2.96, 0.41, 0, 0, 0, 0, 0,
0, 0, 0, 1.42, 0.17, 12.91, 1.22, 42.88, 5.7, 83.86, 9.77, 100,
12.53, 77.91, 10.65, 40.15, 6.08, 10.69, 1.55, 1.47, 0.17, 0,
0, 0, 0, 0, 0, 0.82, 0.09, 8.66, 0.91, 29.55, 3.69, 69.49, 8.98,
100, 12.18, 98.13, 12.04, 66.12, 7.48, 26.81, 3.83, 5.47, 0.64,
1.05, 0.05, 0, 0, 0.46, 0.09, 4.13, 0.7, 20.23, 2.65, 51.74,
5.99, 87.84, 10.86, 100, 13.36, 80.05, 10.81, 44.59, 6.4, 17.31,
1.86, 4.04, 0.6, 0.28, 0.05, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 100, 13.16, 34.55, 3.92, 0.17, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 12.48, 63.63, 8.56, 9.58,
1.47, 0.12, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 100, 14.11, 98.6, 12.13, 31.05, 4.67, 2.93, 0.21, 0.08, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 79.84, 11.22, 100, 12.6,
51.49, 6.67, 11.69, 1.32, 0.85, 0.06, 0.03, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 59.76, 8.6, 100, 13.72, 65.7, 8.43, 19.82,
3.05, 2.86, 0.57, 0.31, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
52.77, 7.04, 100, 12.74, 82.14, 10.06, 34.56, 4.34, 8.21, 0.97,
1.42, 0.16, 0.06, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 44.13, 6.53,
100, 13.32, 98.36, 11.89, 53.27, 6.14, 17.1, 2.5, 3.07, 0.29,
0.39, 0.04, 0, 0, 0, 0, 0, 0, 0, 0, 33.36, 4.74, 90.92, 11.81,
100, 13.16, 64.93, 8.75, 25.77, 3.35, 5.98, 0.77, 0.84, 0.15,
0.04, 0, 0, 0, 0, 0, 0, 0, 0, 0, 26.85, 3.75, 79.08, 10.71, 100,
12.35, 80.84, 10.21, 36.79, 5.62, 10.94, 2.18, 2.52, 0.11, 0.46,
0, 0, 0, 0, 0, 0, 0, 0, 0, 19.82, 2.88, 70.8, 9.23, 100, 14.1,
88.21, 11.51, 48.22, 6.59, 19.66, 2.52, 4, 0.72, 0.96, 0.16,
0.12, 0.04, 0, 0, 0, 0, 0, 0, 0, 0, 76.28, 10.04, 100, 13, 25.11,
2.98, 0.2, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 62.46, 8.28, 100, 12.91, 45.34, 5.93, 6.42, 0.75, 0.1, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 49.91, 6.63, 100, 13.38,
65.58, 8.01, 17.23, 2.19, 1.16, 0.24, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 45.71, 6.38, 100, 14.15, 82.74, 9.31, 32.16, 
4.23, 5.96, 1, 0.51, 0.15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 37.86,
4.85, 98.1, 12.47, 100, 12.31, 49.08, 5.8, 13.56, 1.76, 2.71,
0.31, 0.14, 0.03, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 43.98,
6.08, 100, 12.84, 72.24, 9.07, 13.8, 1.82, 0.1, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 38.4, 5.7, 100, 13.08, 85.68, 11.54,
32, 4.25, 3.9, 0.59, 0.03, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32.33,
4.16, 92.83, 11.75, 100, 13.84, 49.8, 7.07, 13.06, 1.51, 0.75,
0.1, 0.03, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25.09, 3.52, 81.05,
10.67, 100, 12.95, 64.24, 8.32, 22.81, 2.52, 3.62, 0.4, 0.23,
0.03, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 19.96, 2.47, 68.05, 9.28,
100, 12.69, 75.07, 10.88, 33.55, 4.55, 8.62, 0.9, 1.43, 0.1,
0.14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25.5, 3.22, 84.74, 11.17, 100,
12.75, 46.61, 6.27, 7.58, 0.72, 0.03, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 20.61, 2.92, 74.55, 8.85, 100, 13.37, 62, 8.22, 
18.92, 2.35, 1.76, 0.19, 0.03, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 17.07, 2.38, 62.66, 7.03, 100, 13.07, 77.42, 9.14, 28.92,
4.79, 6.5, 0.69, 0.43, 0.03, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 15.18,
2.22, 55.24, 7.26, 100, 13.11, 91.63, 11.66, 48.72, 5.96, 14.22,
2.3, 2.18, 0.3, 0.26, 0, 0, 0, 0, 0, 0, 0, 12.31, 1.55, 51.33,
7.65, 96.93, 12.39, 100, 12.58, 60.41, 8.24, 24.33, 3.5, 5.26,
0.92, 0.76, 0.04, 0.04, 0, 0, 0, 0, 0, 0, 0, 0, 0, 15.67, 2.12,
60.92, 8.03, 100, 13.58, 80.05, 10.89, 32.24, 4.28, 3.45, 0.73,
0.03, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 11.78, 1.79, 54.69, 6.59,
100, 12.11, 94.52, 10.71, 46.24, 5.62, 11.28, 1.61, 1.15, 0.07,
0, 0, 0, 0, 0, 0, 0, 0, 9.69, 1.27, 48.48, 6.18, 95.83, 12.58,
100, 13.93, 61.25, 7.8, 21.11, 3.05, 3.82, 0.54, 0.39, 0.04,
0, 0, 0, 0, 0, 0, 0, 0, 6.76, 1.28, 37.15, 5.72, 83.53, 10.94,
100, 13.18, 71.43, 8.81, 32.86, 3.59, 8.47, 1.04, 1.51, 0.15,
0.08, 0.04, 0, 0, 0, 0, 0, 0, 6.2, 0.58, 28.92, 4.05, 73.89,
9.04, 100, 12.47, 76.73, 10.44, 42.32, 4.33, 14.65, 2.3, 3.12, 
0.19, 0.43, 0.04, 0, 0, 0, 0, 0, 0, 7.15, 1.08, 42.78, 5.17,
88.11, 11.24, 100, 12.9, 60.31, 7.4, 18.46, 2.12, 2.33, 0.14,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8.03, 1, 34.79, 4.65, 82.87, 10.87,
100, 12.67, 78.03, 8.79, 30.07, 3.8, 6.64, 1.08, 0.65, 0.04,
0.04, 0, 0, 0, 0, 0, 0, 0, 4.72, 0.7, 28.35, 3.4, 73.64, 10.58,
100, 13.04, 80.44, 10.89, 43.07, 5.31, 11.68, 1.56, 2.5, 0.43,
0.16, 0, 0, 0, 0, 0, 0, 0, 4.15, 0.59, 23.58, 3.22, 68.29, 9.27,
100, 13.42, 94.79, 12.32, 55.93, 7.49, 21.34, 3.22, 4.06, 0.76,
0.8, 0.13, 0, 0, 0, 0, 2.95, 0.45, 20.2, 2.68, 62.33, 7.15, 97.36,
12.82, 100, 13.99, 68.54, 9.03, 31.68, 3.49, 10.32, 1.47, 2.06,
0.09, 0.13, 0.09, 0, 0), distribution_name = c("Br1", "Br1",
"Br1", "Br1", "Br1", "Br1", "Br1", "Br1", "Br1", "Br1", "Br1",
"Br1", "Br1", "Br1", "Br1", "Br1", "Br1", "Br1", "Br1", "Br1",
"Br1", "Br1", "Br1", "Br1", "Br2", "Br2", "Br2", "Br2", "Br2",
"Br2", "Br2", "Br2", "Br2", "Br2", "Br2", "Br2", "Br2", "Br2",
"Br2", "Br2", "Br2", "Br2", "Br2", "Br2", "Br2", "Br2", "Br2",
"Br2", "Br3", "Br3", "Br3", "Br3", "Br3", "Br3", "Br3", "Br3",
"Br3", "Br3", "Br3", "Br3", "Br3", "Br3", "Br3", "Br3", "Br3", 
"Br3", "Br3", "Br3", "Br3", "Br3", "Br3", "Br3", "Br4", "Br4",
"Br4", "Br4", "Br4", "Br4", "Br4", "Br4", "Br4", "Br4", "Br4",
"Br4", "Br4", "Br4", "Br4", "Br4", "Br4", "Br4", "Br4", "Br4",
"Br4", "Br4", "Br4", "Br4", "Br5", "Br5", "Br5", "Br5", "Br5",
"Br5", "Br5", "Br5", "Br5", "Br5", "Br5", "Br5", "Br5", "Br5",
"Br5", "Br5", "Br5", "Br5", "Br5", "Br5", "Br5", "Br5", "Br5",
"Br5", "Br6", "Br6", "Br6", "Br6", "Br6", "Br6", "Br6", "Br6",
"Br6", "Br6", "Br6", "Br6", "Br6", "Br6", "Br6", "Br6", "Br6",
"Br6", "Br6", "Br6", "Br6", "Br6", "Br6", "Br6", "Br7", "Br7",
"Br7", "Br7", "Br7", "Br7", "Br7", "Br7", "Br7", "Br7", "Br7",
"Br7", "Br7", "Br7", "Br7", "Br7", "Br7", "Br7", "Br7", "Br7",
"Br7", "Br7", "Br7", "Br7", "Br8", "Br8", "Br8", "Br8", "Br8",
"Br8", "Br8", "Br8", "Br8", "Br8", "Br8", "Br8", "Br8", "Br8",
"Br8", "Br8", "Br8", "Br8", "Br8", "Br8", "Br8", "Br8", "Br8",
"Br8", "Br9", "Br9", "Br9", "Br9", "Br9", "Br9", "Br9", "Br9",
"Br9", "Br9", "Br9", "Br9", "Br9", "Br9", "Br9", "Br9", "Br9", 
"Br9", "Br9", "Br9", "Br9", "Br9", "Br9", "Br9", "Br10", "Br10",
"Br10", "Br10", "Br10", "Br10", "Br10", "Br10", "Br10", "Br10",
"Br10", "Br10", "Br10", "Br10", "Br10", "Br10", "Br10", "Br10",
"Br10", "Br10", "Br10", "Br10", "Br10", "Br10", "Cl1", "Cl1",
"Cl1", "Cl1", "Cl1", "Cl1", "Cl1", "Cl1", "Cl1", "Cl1", "Cl1",
"Cl1", "Cl1", "Cl1", "Cl1", "Cl1", "Cl1", "Cl1", "Cl1", "Cl1",
"Cl1", "Cl1", "Cl1", "Cl1", "Cl2", "Cl2", "Cl2", "Cl2", "Cl2",
"Cl2", "Cl2", "Cl2", "Cl2", "Cl2", "Cl2", "Cl2", "Cl2", "Cl2",
"Cl2", "Cl2", "Cl2", "Cl2", "Cl2", "Cl2", "Cl2", "Cl2", "Cl2",
"Cl2", "Cl3", "Cl3", "Cl3", "Cl3", "Cl3", "Cl3", "Cl3", "Cl3",
"Cl3", "Cl3", "Cl3", "Cl3", "Cl3", "Cl3", "Cl3", "Cl3", "Cl3",
"Cl3", "Cl3", "Cl3", "Cl3", "Cl3", "Cl3", "Cl3", "Cl4", "Cl4",
"Cl4", "Cl4", "Cl4", "Cl4", "Cl4", "Cl4", "Cl4", "Cl4", "Cl4",
"Cl4", "Cl4", "Cl4", "Cl4", "Cl4", "Cl4", "Cl4", "Cl4", "Cl4",
"Cl4", "Cl4", "Cl4", "Cl4", "Cl5", "Cl5", "Cl5", "Cl5", "Cl5",
"Cl5", "Cl5", "Cl5", "Cl5", "Cl5", "Cl5", "Cl5", "Cl5", "Cl5",
"Cl5", "Cl5", "Cl5", "Cl5", "Cl5", "Cl5", "Cl5", "Cl5", "Cl5", 
"Cl5", "Cl6", "Cl6", "Cl6", "Cl6", "Cl6", "Cl6", "Cl6", "Cl6",
"Cl6", "Cl6", "Cl6", "Cl6", "Cl6", "Cl6", "Cl6", "Cl6", "Cl6",
"Cl6", "Cl6", "Cl6", "Cl6", "Cl6", "Cl6", "Cl6", "Cl7", "Cl7",
"Cl7", "Cl7", "Cl7", "Cl7", "Cl7", "Cl7", "Cl7", "Cl7", "Cl7",
"Cl7", "Cl7", "Cl7", "Cl7", "Cl7", "Cl7", "Cl7", "Cl7", "Cl7",
"Cl7", "Cl7", "Cl7", "Cl7", "Cl8", "Cl8", "Cl8", "Cl8", "Cl8",
"Cl8", "Cl8", "Cl8", "Cl8", "Cl8", "Cl8", "Cl8", "Cl8", "Cl8",
"Cl8", "Cl8", "Cl8", "Cl8", "Cl8", "Cl8", "Cl8", "Cl8", "Cl8",
"Cl8", "Cl9", "Cl9", "Cl9", "Cl9", "Cl9", "Cl9", "Cl9", "Cl9",
"Cl9", "Cl9", "Cl9", "Cl9", "Cl9", "Cl9", "Cl9", "Cl9", "Cl9",
"Cl9", "Cl9", "Cl9", "Cl9", "Cl9", "Cl9", "Cl9", "Cl10", "Cl10",
"Cl10", "Cl10", "Cl10", "Cl10", "Cl10", "Cl10", "Cl10", "Cl10",
"Cl10", "Cl10", "Cl10", "Cl10", "Cl10", "Cl10", "Cl10", "Cl10", 
"Cl10", "Cl10", "Cl10", "Cl10", "Cl10", "Cl10", "Br1Cl1", "Br1Cl1",
"Br1Cl1", "Br1Cl1", "Br1Cl1", "Br1Cl1", "Br1Cl1", "Br1Cl1", "Br1Cl1",
"Br1Cl1", "Br1Cl1", "Br1Cl1", "Br1Cl1", "Br1Cl1", "Br1Cl1", "Br1Cl1",
"Br1Cl1", "Br1Cl1", "Br1Cl1", "Br1Cl1", "Br1Cl1", "Br1Cl1", "Br1Cl1",
"Br1Cl1", "Br1Cl2", "Br1Cl2", "Br1Cl2", "Br1Cl2", "Br1Cl2", "Br1Cl2",
"Br1Cl2", "Br1Cl2", "Br1Cl2", "Br1Cl2", "Br1Cl2", "Br1Cl2", "Br1Cl2",
"Br1Cl2", "Br1Cl2", "Br1Cl2", "Br1Cl2", "Br1Cl2", "Br1Cl2", "Br1Cl2",
"Br1Cl2", "Br1Cl2", "Br1Cl2", "Br1Cl2", "Br1Cl3", "Br1Cl3", "Br1Cl3",
"Br1Cl3", "Br1Cl3", "Br1Cl3", "Br1Cl3", "Br1Cl3", "Br1Cl3", "Br1Cl3",
"Br1Cl3", "Br1Cl3", "Br1Cl3", "Br1Cl3", "Br1Cl3", "Br1Cl3", "Br1Cl3",
"Br1Cl3", "Br1Cl3", "Br1Cl3", "Br1Cl3", "Br1Cl3", "Br1Cl3", "Br1Cl3",
"Br1Cl4", "Br1Cl4", "Br1Cl4", "Br1Cl4", "Br1Cl4", "Br1Cl4", "Br1Cl4",
"Br1Cl4", "Br1Cl4", "Br1Cl4", "Br1Cl4", "Br1Cl4", "Br1Cl4", "Br1Cl4",
"Br1Cl4", "Br1Cl4", "Br1Cl4", "Br1Cl4", "Br1Cl4", "Br1Cl4", "Br1Cl4",
"Br1Cl4", "Br1Cl4", "Br1Cl4", "Br1Cl5", "Br1Cl5", "Br1Cl5", "Br1Cl5",
"Br1Cl5", "Br1Cl5", "Br1Cl5", "Br1Cl5", "Br1Cl5", "Br1Cl5", "Br1Cl5",
"Br1Cl5", "Br1Cl5", "Br1Cl5", "Br1Cl5", "Br1Cl5", "Br1Cl5", "Br1Cl5", 
"Br1Cl5", "Br1Cl5", "Br1Cl5", "Br1Cl5", "Br1Cl5", "Br1Cl5", "Br2Cl1",
"Br2Cl1", "Br2Cl1", "Br2Cl1", "Br2Cl1", "Br2Cl1", "Br2Cl1", "Br2Cl1",
"Br2Cl1", "Br2Cl1", "Br2Cl1", "Br2Cl1", "Br2Cl1", "Br2Cl1", "Br2Cl1",
"Br2Cl1", "Br2Cl1", "Br2Cl1", "Br2Cl1", "Br2Cl1", "Br2Cl1", "Br2Cl1",
"Br2Cl1", "Br2Cl1", "Br2Cl2", "Br2Cl2", "Br2Cl2", "Br2Cl2", "Br2Cl2",
"Br2Cl2", "Br2Cl2", "Br2Cl2", "Br2Cl2", "Br2Cl2", "Br2Cl2", "Br2Cl2",
"Br2Cl2", "Br2Cl2", "Br2Cl2", "Br2Cl2", "Br2Cl2", "Br2Cl2", "Br2Cl2",
"Br2Cl2", "Br2Cl2", "Br2Cl2", "Br2Cl2", "Br2Cl2", "Br2Cl3", "Br2Cl3",
"Br2Cl3", "Br2Cl3", "Br2Cl3", "Br2Cl3", "Br2Cl3", "Br2Cl3", "Br2Cl3",
"Br2Cl3", "Br2Cl3", "Br2Cl3", "Br2Cl3", "Br2Cl3", "Br2Cl3", "Br2Cl3",
"Br2Cl3", "Br2Cl3", "Br2Cl3", "Br2Cl3", "Br2Cl3", "Br2Cl3", "Br2Cl3",
"Br2Cl3", "Br2Cl4", "Br2Cl4", "Br2Cl4", "Br2Cl4", "Br2Cl4", "Br2Cl4",
"Br2Cl4", "Br2Cl4", "Br2Cl4", "Br2Cl4", "Br2Cl4", "Br2Cl4", "Br2Cl4",
"Br2Cl4", "Br2Cl4", "Br2Cl4", "Br2Cl4", "Br2Cl4", "Br2Cl4", "Br2Cl4",
"Br2Cl4", "Br2Cl4", "Br2Cl4", "Br2Cl4", "Br2Cl5", "Br2Cl5", "Br2Cl5", 
"Br2Cl5", "Br2Cl5", "Br2Cl5", "Br2Cl5", "Br2Cl5", "Br2Cl5", "Br2Cl5",
"Br2Cl5", "Br2Cl5", "Br2Cl5", "Br2Cl5", "Br2Cl5", "Br2Cl5", "Br2Cl5",
"Br2Cl5", "Br2Cl5", "Br2Cl5", "Br2Cl5", "Br2Cl5", "Br2Cl5", "Br2Cl5",
"Br3Cl1", "Br3Cl1", "Br3Cl1", "Br3Cl1", "Br3Cl1", "Br3Cl1", "Br3Cl1",
"Br3Cl1", "Br3Cl1", "Br3Cl1", "Br3Cl1", "Br3Cl1", "Br3Cl1", "Br3Cl1",
"Br3Cl1", "Br3Cl1", "Br3Cl1", "Br3Cl1", "Br3Cl1", "Br3Cl1", "Br3Cl1",
"Br3Cl1", "Br3Cl1", "Br3Cl1", "Br3Cl2", "Br3Cl2", "Br3Cl2", "Br3Cl2",
"Br3Cl2", "Br3Cl2", "Br3Cl2", "Br3Cl2", "Br3Cl2", "Br3Cl2", "Br3Cl2",
"Br3Cl2", "Br3Cl2", "Br3Cl2", "Br3Cl2", "Br3Cl2", "Br3Cl2", "Br3Cl2",
"Br3Cl2", "Br3Cl2", "Br3Cl2", "Br3Cl2", "Br3Cl2", "Br3Cl2", "Br3Cl3",
"Br3Cl3", "Br3Cl3", "Br3Cl3", "Br3Cl3", "Br3Cl3", "Br3Cl3", "Br3Cl3",
"Br3Cl3", "Br3Cl3", "Br3Cl3", "Br3Cl3", "Br3Cl3", "Br3Cl3", "Br3Cl3",
"Br3Cl3", "Br3Cl3", "Br3Cl3", "Br3Cl3", "Br3Cl3", "Br3Cl3", "Br3Cl3",
"Br3Cl3", "Br3Cl3", "Br3Cl4", "Br3Cl4", "Br3Cl4", "Br3Cl4", "Br3Cl4", 
"Br3Cl4", "Br3Cl4", "Br3Cl4", "Br3Cl4", "Br3Cl4", "Br3Cl4", "Br3Cl4",
"Br3Cl4", "Br3Cl4", "Br3Cl4", "Br3Cl4", "Br3Cl4", "Br3Cl4", "Br3Cl4",
"Br3Cl4", "Br3Cl4", "Br3Cl4", "Br3Cl4", "Br3Cl4", "Br3Cl5", "Br3Cl5",
"Br3Cl5", "Br3Cl5", "Br3Cl5", "Br3Cl5", "Br3Cl5", "Br3Cl5", "Br3Cl5",
"Br3Cl5", "Br3Cl5", "Br3Cl5", "Br3Cl5", "Br3Cl5", "Br3Cl5", "Br3Cl5",
"Br3Cl5", "Br3Cl5", "Br3Cl5", "Br3Cl5", "Br3Cl5", "Br3Cl5", "Br3Cl5",
"Br3Cl5", "Br4Cl1", "Br4Cl1", "Br4Cl1", "Br4Cl1", "Br4Cl1", "Br4Cl1",
"Br4Cl1", "Br4Cl1", "Br4Cl1", "Br4Cl1", "Br4Cl1", "Br4Cl1", "Br4Cl1",
"Br4Cl1", "Br4Cl1", "Br4Cl1", "Br4Cl1", "Br4Cl1", "Br4Cl1", "Br4Cl1",
"Br4Cl1", "Br4Cl1", "Br4Cl1", "Br4Cl1", "Br4Cl2", "Br4Cl2", "Br4Cl2",
"Br4Cl2", "Br4Cl2", "Br4Cl2", "Br4Cl2", "Br4Cl2", "Br4Cl2", "Br4Cl2",
"Br4Cl2", "Br4Cl2", "Br4Cl2", "Br4Cl2", "Br4Cl2", "Br4Cl2", "Br4Cl2",
"Br4Cl2", "Br4Cl2", "Br4Cl2", "Br4Cl2", "Br4Cl2", "Br4Cl2", "Br4Cl2",
"Br4Cl3", "Br4Cl3", "Br4Cl3", "Br4Cl3", "Br4Cl3", "Br4Cl3", "Br4Cl3", 
"Br4Cl3", "Br4Cl3", "Br4Cl3", "Br4Cl3", "Br4Cl3", "Br4Cl3", "Br4Cl3",
"Br4Cl3", "Br4Cl3", "Br4Cl3", "Br4Cl3", "Br4Cl3", "Br4Cl3", "Br4Cl3",
"Br4Cl3", "Br4Cl3", "Br4Cl3", "Br4Cl4", "Br4Cl4", "Br4Cl4", "Br4Cl4",
"Br4Cl4", "Br4Cl4", "Br4Cl4", "Br4Cl4", "Br4Cl4", "Br4Cl4", "Br4Cl4",
"Br4Cl4", "Br4Cl4", "Br4Cl4", "Br4Cl4", "Br4Cl4", "Br4Cl4", "Br4Cl4",
"Br4Cl4", "Br4Cl4", "Br4Cl4", "Br4Cl4", "Br4Cl4", "Br4Cl4", "Br4Cl5",
"Br4Cl5", "Br4Cl5", "Br4Cl5", "Br4Cl5", "Br4Cl5", "Br4Cl5", "Br4Cl5",
"Br4Cl5", "Br4Cl5", "Br4Cl5", "Br4Cl5", "Br4Cl5", "Br4Cl5", "Br4Cl5",
"Br4Cl5", "Br4Cl5", "Br4Cl5", "Br4Cl5", "Br4Cl5", "Br4Cl5", "Br4Cl5",
"Br4Cl5", "Br4Cl5", "Br5Cl1", "Br5Cl1", "Br5Cl1", "Br5Cl1", "Br5Cl1",
"Br5Cl1", "Br5Cl1", "Br5Cl1", "Br5Cl1", "Br5Cl1", "Br5Cl1", "Br5Cl1",
"Br5Cl1", "Br5Cl1", "Br5Cl1", "Br5Cl1", "Br5Cl1", "Br5Cl1", "Br5Cl1",
"Br5Cl1", "Br5Cl1", "Br5Cl1", "Br5Cl1", "Br5Cl1", "Br5Cl2", "Br5Cl2",
"Br5Cl2", "Br5Cl2", "Br5Cl2", "Br5Cl2", "Br5Cl2", "Br5Cl2", "Br5Cl2",
"Br5Cl2", "Br5Cl2", "Br5Cl2", "Br5Cl2", "Br5Cl2", "Br5Cl2", "Br5Cl2", 
"Br5Cl2", "Br5Cl2", "Br5Cl2", "Br5Cl2", "Br5Cl2", "Br5Cl2", "Br5Cl2",
"Br5Cl2", "Br5Cl3", "Br5Cl3", "Br5Cl3", "Br5Cl3", "Br5Cl3", "Br5Cl3",
"Br5Cl3", "Br5Cl3", "Br5Cl3", "Br5Cl3", "Br5Cl3", "Br5Cl3", "Br5Cl3",
"Br5Cl3", "Br5Cl3", "Br5Cl3", "Br5Cl3", "Br5Cl3", "Br5Cl3", "Br5Cl3",
"Br5Cl3", "Br5Cl3", "Br5Cl3", "Br5Cl3", "Br5Cl4", "Br5Cl4", "Br5Cl4",
"Br5Cl4", "Br5Cl4", "Br5Cl4", "Br5Cl4", "Br5Cl4", "Br5Cl4", "Br5Cl4",
"Br5Cl4", "Br5Cl4", "Br5Cl4", "Br5Cl4", "Br5Cl4", "Br5Cl4", "Br5Cl4",
"Br5Cl4", "Br5Cl4", "Br5Cl4", "Br5Cl4", "Br5Cl4", "Br5Cl4", "Br5Cl4",
"Br5Cl5", "Br5Cl5", "Br5Cl5", "Br5Cl5", "Br5Cl5", "Br5Cl5", "Br5Cl5",
"Br5Cl5", "Br5Cl5", "Br5Cl5", "Br5Cl5", "Br5Cl5", "Br5Cl5", "Br5Cl5",
"Br5Cl5", "Br5Cl5", "Br5Cl5", "Br5Cl5", "Br5Cl5", "Br5Cl5", "Br5Cl5",
"Br5Cl5", "Br5Cl5", "Br5Cl5")), class = "data.frame", row.names = c(NA,
-1080L))

### Split ChromaTOF-exported MSP file into list of individual spectra

# Terms: 'sample' is a single MSP file generated from the analysis of a sample. Process all sample MSP files in the folder, where input_sample_list is a vector of paths to the input MSP spectra files.
input_sample_list <- list.files("Input Spectra", full.names = TRUE) 

# for (i in 1:length(input_sample_list)) { # Increment through path list
# TODO complete all tasks for single file then use rm(list()) to clear memory except for functions before going to the next file?
i <- 1 # For dev: work with single file
working_sample <- input_sample_list[i] # Select path from list
tmp_name <- strsplit(working_sample, split = '/', fixed = TRUE)[[1]][2] # Select file name from path
working_sample_name <- substr(tmp_name, 1, nchar(tmp_name) - 4) # Remove '.txt' from file name # TODO Use or delete?
working_sample_spectra <- readLines(working_sample) # Read file in to character vector 

cat('Processing ', tmp_name, '... \n', sep = '')

# Make a vector (spectrum_factor) that defines the separate spectra within the MSP file, where spectrum_factor is defined by spectrum_count, which iterates by + 1 if it encounters a newline that defines the next spectrum. 
spectrum_factor <- vector(mode = "integer", length = length(working_sample_spectra))
spectrum_count <- 1
for(i in 1: length(working_sample_spectra)) {
  if(working_sample_spectra[i] != "") { 
    spectrum_factor[i] <- spectrum_count 
  } else {
    spectrum_factor[i] <- spectrum_count
    spectrum_count <- spectrum_count + 1 
  }
}

# Using spectrum_factor, split the spectra into list elements as defined by spectrum_factor and convert to a data frame.
spectrum_list <- split(working_sample_spectra, f = spectrum_factor)

# TODO If end of file contains more than one newline, ignore
# REPORT number of spectra in file

# Preallocate the data frame size (N) with an estimate of the total number of m/z peaks in the file: find the number of rows in each spectrum, exclude 3 rows per spectrum for the headers and blank line, find the sum total, and multiply by 10 (the max number of allowed m/z peaks per row). Create the empty data frame (df).
N <- sum(lengths(spectrum_list) - 3) * 10 
df <- data.frame(spectrum_number = rep(NA, N), subnominal_mz = rep(NA, N), nominal_mz = rep(NA, N), intensity = rep(NA, N))

# Track data frame row number (row_num) and iterate through the individual spectra (spectrum_list_i) from spectrum_list. Split out peaks (m/z and intensity pairs) into list subelements (x). Remove `character(0)` elements from list that results from newlines.
n <- length(spectrum_list)
row_num <- 1 
cat('Total number of spectra:', n, '\n')
cat('Reading in spectrum:\n')

for(spectrum_list_i in 1:length(spectrum_list)) { 

  cat(sprintf('\r'), spectrum_list_i, sep = '')

  x <- strsplit(spectrum_list[[spectrum_list_i]], split = ";", fixed = TRUE)
  x <- x[lengths(x) != 0] 

  # Within an individual spectrum, iterate through the rows (i), skiping NAME and Num Peaks. Then iterate through the m/z peaks (m/z and intensity pairs) within a row (j) and make a vector of the pair (mz_int). Remove if empty and use round() to deterine the nominal m/z. Construct a list (a) of the values and add this as a row to df.
  for(i in 3:length(x)) { 
    for(j in 1:length(x[[i]])) { 
      mz_int <- strsplit(x[[i]][j], split = "[[:space:]]")[[1]] 
      mz_int <- mz_int[mz_int != ""] 
      nominal_mz <- round(as.numeric(mz_int[1]))
      a <- list(as.numeric(spectrum_list_i), as.numeric(mz_int)[1], as.numeric(nominal_mz), as.numeric(mz_int)[2]) 
      df[row_num, ] <- a
      row_num <- row_num + 1

    }
  }
}

### Check for duplicated nominal m/z values and combine

# Remove unused rows (NA). For the remaining data, sum the intensities of duplicated nominal_mz values within the spectrum. The subnominal_mz for the the m/z peak with the larger original intensity is kept, but with the original intensity replaced by the summed intensity. The smaller intensity m/z peak is deleted. Columns are then reordered.
df <- na.omit(df) 
df_aggregated <- aggregate(df['intensity'], by = df[c('spectrum_number', 'nominal_mz')], sum)
df_combined <- merge(df, df_aggregated, by = c('spectrum_number', 'nominal_mz'))
df_ordered <- df_combined[order(df_combined$spectrum_number, df_combined$nominal_mz, -df_combined$intensity.x), ] 
df_read <- df_ordered[!duplicated(df_ordered[c('spectrum_number', 'nominal_mz', 'intensity.y')]), ]
rownames(df_read) <- NULL
df_read <- df_read[, c('spectrum_number', 'subnominal_mz', 'nominal_mz', 'intensity.y')]
names(df_read)[names(df_read) == 'intensity.y'] <- 'intensity'

# TODO Make Intensity threshold settable
# REPORT status to separate output file and terminal. 
# REPORT number of dulicated nominal m/z peaks per spectrum
# REPORT number of m/z peaks under intensity threshold
# TODO Double check with more complex examples

### Remove noise spectra

rm(list = setdiff(ls(), c("df_read", "t_distributions")))
cat(' -> Completed...\n')
cat('Memory used to read data frame:', format(object.size(df_read), units = 'Mb'), '\n')

# Remove spectra with less than 10 ions above m/z 100. These are 'noise spectra'. 
n_removed <- 0
for(i in unique(df_read$spectrum_number)) {
  tmp <- df_read[df_read$spectrum_number == i, ]
  if(nrow(tmp[tmp$nominal_mz > 100, ]) < 10) {
    df_read <- df_read[df_read$spectrum_number != i, ]
    n_removed <- n_removed + 1
  }    
}
cat('Number of noise spectra removed:', n_removed, '\n')
cat('New total number of spectra:', length(unique(df_read$spectrum_number)), '\n')

### Make data frame for input into filter set

# Add in placeholders for all nominal m/z values
# Construct df with all nominal mz values and merge. Find min and max values for each spectrum_number and build df using rep().  
cat('Adding in missing nominal m/z values... \n')
z <- data.frame(NULL)
for(i in unique(df_read$spectrum_number)) {
  tmp <- df_read[df_read$spectrum_number == i, ]
  x_min <- min(tmp$nominal_mz)
  x_max <- max(tmp$nominal_mz)
  nominal_values <- x_min:x_max
  y <- data.frame(spectrum_number = i, nominal_mz = nominal_values)
  z <- rbind(z, y)
}

df <- merge(df_read, z, by = c('spectrum_number', 'nominal_mz'), all = TRUE)
df$intensity[is.na(df$intensity)] <- 0
rownames(df_read) <- NULL
cat('Memory used by full nominal m/z data frame:', format(object.size(df), units = 'Mb'), '\n')

#write.csv(df_read, file = 'Output.txt')

### Halogenation filter 1
# TODO Cite paper
rm(list = setdiff(ls(), c("df", "t_distributions")))

# Normalize peak intensity to percentage of max peak and only examine m/z > 100 amu.
cat('Applying filter 1 to spectrum (original ID# is printed):\n')
df$filter_1 <- FALSE
y <- data.frame(NULL)

for(i in unique(df$spectrum_number)) {
  cat(sprintf('\r'), i, sep = '')
  x <- df[df$spectrum_number == i, ]
  x$intensity <- with(x, intensity / max(intensity) * 100)
  x <- x[x$nominal_mz > 100, ]

  for(i in (nrow(x) - 6):5) {
    
    if(x$intensity[i] > 3) {
      
      # Calculate peak ratios.
      
      peakRatio1 <- x$intensity[i + 1] / x$intensity[i]
      peakRatio2 <- x$intensity[i + 2] / x$intensity[i]
      peakRatio3 <- x$intensity[i + 3] / x$intensity[i]
      peakRatio4 <- x$intensity[i + 4] / x$intensity[i]
      peakRatio5 <- x$intensity[i + 5] / x$intensity[i]
      peakRatio6 <- x$intensity[i + 6] / x$intensity[i]
      peakRatioa <- x$intensity[i - 1] / x$intensity[i]
      peakRatiob <- x$intensity[i - 2] / x$intensity[i]
      peakRatioc <- x$intensity[i - 4] / x$intensity[i]
      
      # Determine if peak ratios pass rules.
      
      if(peakRatio1 < 0.5 & peakRatio1 > 0.003 &   # Rule 1
         peakRatio2 > 0.225 & peakRatio2 <= 1 &    # Rule 2
         peakRatio1 < peakRatio2 &                 # Rule 3
         peakRatio3 < 0.5 & peakRatio3 > 0.001 &   # Rule 4
         peakRatio3 < peakRatio1 &                 # Rule 5
         peakRatio4 < peakRatio2 &                 # Rule 6
         peakRatio5 < peakRatio1 &                 # Rule 7
         peakRatio6 < peakRatio2 &                 # Rule 8
         peakRatioa < 0.5 &                        # Rule 9
         peakRatiob < 1 &                          # Rule 10
         peakRatioc < 1) {                         # Rule 11
        
        x$filter_1[i] <- TRUE 
      }
    } 
  } 
  x <- x[, c('spectrum_number', 'nominal_mz', 'filter_1')]
  y <- rbind(y, x)

}


df <- merge(df, y, by = c('spectrum_number', 'nominal_mz')) # df holds the results of filter 1, select columns
df <- subset(df, select = -filter_1.x) 
names(df)[names(df) == 'filter_1.y'] <- 'filter_1' # rename column from merge
df <- df[order(df$spectrum_number, df$nominal_mz), ]
cat(' -> Completed...\n') 
rm(list = setdiff(ls(), c("df", "t_distributions")))
cat('Total number of halogenated spectra per filter 1:', length(unique(df$spectrum_number[df$filter_1 == TRUE])), 'out of', length(unique(df$spectrum_number)), '\n')
cat('Memory used by data frame with filter 1 results:', format(object.size(df), units = 'Mb'), '\n')

### Halogenation filter 2

# Iterate through the passing filter 1 results. For each, determine the m/z range for selection and extract the experimental distribution. Convert intensities to the percent intensity within the distribution. This conversion also changes intensity values that are NA to percent_intensity values that are 0.

filter_1 <- df[df$filter_1 == TRUE,] # passing rows of filter 1
distribution_names <- unique(t_distributions$distribution_name) # theoretical distribution names
filter_2 <- data.frame(NULL) # to collect filter 2 results
monoisotopic_peaks <- data.frame(NULL) # to collect monoisotopic peaks

cat('Applying filter 2 to spectrum (count, original ID#):\n')
m <- nrow(filter_1)

for(i in 1:nrow(filter_1)) {

    n <- filter_1$spectrum_number[i]
    cat(sprintf('\r'), i, " of ", m, ", ", n, sep = '')

    mz_range <- data.frame(spectrum_number = n, nominal_mz = seq(from = filter_1$nominal_mz[i] - 10, to = filter_1$nominal_mz[i] + 13, by = 1)) # set up mz range
    e_distribution <- merge(mz_range, df, by = c('spectrum_number', 'nominal_mz'), all.x = TRUE, all.y = FALSE) # fill in missing rows in the experimental distribution
    e_distribution$intensity[is.na(e_distribution$intensity)] <- 0

    # Calculate a similarity score for the experimental distribution against every theoretical distribution. Both have 24 nominal m/z values with 10 below the maximum intensity peak and 13 above.

    match_results <- data.frame(NULL)

    for(j in 1:length(distribution_names)) {

      t_distribution <- t_distributions[t_distributions$distribution_name == distribution_names[j], ] # theoretical distribution
      match_df <- cbind(e_distribution, t_distribution) # bound distributions
      match_df <- match_df[!is.na(match_df$theoretical_mz), ] # exclude m/z range outside of theoretical distribution
      match_df$percent_intensity <- with(match_df, intensity / max(intensity) * 100) # scale the experimental distribution intensities
      u <- match_df$percent_intensity # experimental intensity vector
      v <- match_df$theoretical_percent_intensity # theoretical intensity vector
      similarity_score <- as.vector((u %*% v) / (sqrt(sum(u^2)) * sqrt(sum(v^2))))
      match_df$similarity_score <- similarity_score # match results for a single theoretical distribution
      match_df$monoisotopic_peak <- c(TRUE, rep(FALSE, nrow(match_df) - 1)) # label the monoisotopic peak
      match_results <- rbind(match_results, match_df) # collect all similarity score results

    }

    best_match <- match_results[match_results$similarity_score == max(match_results$similarity_score), ] # select the best match
    filter_2 <- rbind(filter_2, best_match) # collect the best match records in the sample

}

filter_2 <- filter_2[c('spectrum_number', 'nominal_mz', 'percent_intensity', 'theoretical_percent_intensity', 'distribution_name', 'similarity_score', 'monoisotopic_peak')]
cat(' -> Completed...\n') 

df <- merge(df, filter_2, by = c('spectrum_number', 'nominal_mz'), all.x = TRUE, all.y = FALSE)
cat('Total number of halogenated spectra passing filter 2:', length(unique(filter_2$spectrum_number)), '\n') # TODO apply threshold, currently the same as for filter 1
cat('Memory used by data frame with filter 1 and 2 results:', format(object.size(df), units = 'Mb'), '\n')


for k in {01..24};do
	j="LG"$k
	head -n 1 HP.annotation.named.gff > HP.annotation.named.$j.gff
	grep ^$j HP.annotation.named.gff >> HP.annotation.named.$j.gff
	
	gzip HP.annotation.named.$j.gff
done

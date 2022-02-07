# Procesamiento reads MinION

# Consideraciones:
- El contenido entre corchetes se modificará en cada análisis, corresponde a los nombres de las muestras.
- Crear una carpeta donde se va a llamar al paso #1.
- La carpeta principal va a contener una carpeta por cada muestra a procesar. Reproducir los pasos adentro del directorio de cada una de las muestras.

# Workflow:

### Preprocesamiento de datos:

1. Ingresar al terminar de Ubuntu.

2. Activar el entorno Conda:
```
conda activate artic-ncov2019
```
3. Acceder a la carpeta de artic:
```
cd ~/artic-ncov2019
```
4. Seleccionar los archivos fastq_pass y enviarlos a Rampart
```
rampart --protocol rampart/ --basecalledPath ~/../genomica_inspi/../../mnt/e/[carpeta_de_interes]/fastq_pass ––clearAnnotated
```
5.	Abrir en el navegador localhost:3000 para ver barcodes

6.	Exportar cada barcode de forma manual

7.	Los archivos se guardarán en cd ~/artic-ncov2019

8.	Generar una carpeta en Windows llamada files_rampart_ready

9.	Copiar los archivos guardados en cd ~/artic-ncov2019 con extensión .fastq

```
mv *.fastq /mnt/e/[carpeta_de_interes]/files_rampart_ready
```
10.	En la carpeta files_rampart_ready generar carpetas para cada barcode. Cada carpeta contendrá 3 archivos: el binned_barcode con extensión .fastq que salió de RAMPART, el archivo references.fa (genoma de referencia) y archivo Python renameHeader.

### Detección de variantes:
1. Indexar el genoma de referencia
```
./bwa index [genomaRef.fa] 
  samtools faidx [genomaRef.fa]
```

2. Alineamiento
```
./bwa mem [genomaRef.fa][muestraBarcode.fa] > [output.sam]
```

3. Conversión de formato SAM a BAM
```
samtools view -S -b [output.sam] > [output.bam]
```

4. Ordenar los alineamientos
```
samtools sort -o [output_sorted.bam] [output.bam] 
```

5. Indexar el archivo BAM
```
samtools index [output_sorted.bam]
```

6. Calcular la cobertura por base
```
mosdepth [prefijo] [output_sorted.bam]
```

7. Extraer bases con 0 cobertura
```
zcat [prefijo.per-base.bed.gz] | awk '$4==0 {print}' > zero-coverage.bed
```

8.  Marcar las bases que no tuvieron cobertura
```
bedtools maskfasta -fi [references.fa] -bed [zero-coverage.bed] -fo [maskedReference.fa]
```

9. Calcular la cobertura de las lecturas
```
bcftools mpileup -O b -o [raw.bcf] -f [maskedReference.fa] [output_sorted.bam]
```

10. Detectar el polimorfismo de nucleótido único
```
bcftools call --ploidy 1 -m -v -o [variant.vcf] [raw.bcf]
```

11. Comprimir [variant.vcf] y crear su archivo indexado
```
bgzip -c [variant.vcf] > [variant.vcf.gz]
tabix -p vcf [variant.vcf.gz]
```

12. Crear la secuencia de consenso
```
cat [maskedReference.fa] | bcftools consensus [variant.vcf.gz] > consensus.fa
```

13. Cambiar nombre del encabezado y archivo
```
renameHeader.py [consensus.fa] [nombre archivo y encabezado]
```

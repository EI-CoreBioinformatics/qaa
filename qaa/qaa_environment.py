from os.path import join



class QAA_Environment(object):
	def __init__(self, config):		
		self.cwd = config.get("cwd", ".")

		# output_dir = Analysis-dir
		self.output_dir = config.get("out_dir", ".")
		# ensure to use absolute paths (I think because of busco?)
		if not self.output_dir.startswith("/"):
			self.output_dir = join(self.cwd, self.output_dir)

		self.qa_dir = self.qaa_dir = join(self.output_dir, "qaa", "survey" if config["survey_assembly"] else "asm")
		self.qc_dir = join(self.output_dir, "qc")
		self.log_dir = join(self.qaa_dir, "log")

		# busco environment
		self.busco_dir = join(self.qaa_dir, "busco")
		self.busco_geno_dir = join(self.busco_dir, "geno")
		self.busco_tran_dir = join(self.busco_dir, "tran")
		self.busco_prot_dir = join(self.busco_dir, "prot")
		self.busco_init_dir = join(config["etc"], "util", "busco_init_dir")
		self.busco_data_dir = config["resources"]["busco_databases"]
		
		# tool dirs
		self.quast_dir = join(self.qaa_dir, "quast")
		self.blob_dir = join(self.qaa_dir, "blobtools")
		self.qualimap_dir = join(self.qaa_dir, "qualimap")
		self.blast_dir = join(self.qaa_dir, "blast")
		self.bam_dir = join(self.qaa_dir, "bam")
		self.prokka_dir = join(self.output_dir, "annotation", "prokka")

package uk.ac.ebi.ena.sra.cram.format.text;

import java.util.List;

import uk.ac.ebi.ena.sra.cram.format.ReadFeature;

 interface ReadFeaturesFormat {

	public void addFeaturesToStringBuilder(Iterable<ReadFeature> features, int readLength,
			StringBuilder sb);

	public List<ReadFeature> asReadFeatureList(String string);
}

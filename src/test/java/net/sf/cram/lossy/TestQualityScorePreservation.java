package net.sf.cram.lossy;

import static org.hamcrest.core.Is.*;
import static org.junit.Assert.*;

import java.util.List;

import org.junit.Test;

public class TestQualityScorePreservation {

	@Test
	public void test1() {
		QualityScorePreservation p = new QualityScorePreservation("m999_8");
		List<PreservationPolicy> policies = p.getPreservationPolicies();

		assertNotNull(p);
		assertEquals(policies.size(), 1);

		PreservationPolicy policy0 = policies.get(0);
		assertThat(policy0.readCategory.type,
				is(ReadCategoryType.LOWER_MAPPING_SCORE));

		assertThat(policy0.readCategory.param, is(999));

		if (policy0.baseCategories != null)
			assertEquals(policy0.baseCategories.isEmpty(), true);

		QualityScoreTreatment treatment = policy0.treatment;
		assertNotNull(treatment);

		assertThat(treatment.type, is(QualityScoreTreatmentType.BIN));
		assertThat(treatment.param, is(8));
	}

	@Test
	public void test2() {
		QualityScorePreservation p = new QualityScorePreservation("R8-N40");
		List<PreservationPolicy> policies = p.getPreservationPolicies();

		assertNotNull(p);
		assertEquals(policies.size(), 2);

		{
			PreservationPolicy policy0 = policies.get(0);
			assertNull(policy0.readCategory);

			List<BaseCategory> baseCategories = policy0.baseCategories;
			assertNotNull(baseCategories);
			assertEquals(baseCategories.size(), 1);

			BaseCategory c0 = baseCategories.get(0);
			assertEquals(c0.type, BaseCategoryType.MATCH);
			assertEquals(c0.param, -1);

			QualityScoreTreatment treatment = policy0.treatment;
			assertNotNull(treatment);

			assertThat(treatment.type, is(QualityScoreTreatmentType.BIN));
			assertThat(treatment.param, is(8));
		}

		{
			PreservationPolicy policy1 = policies.get(1);
			assertNull(policy1.readCategory);

			List<BaseCategory> baseCategories = policy1.baseCategories;
			assertNotNull(baseCategories);
			assertEquals(baseCategories.size(), 1);
			
			BaseCategory c0 = baseCategories.get(0);
			assertEquals(c0.type, BaseCategoryType.MISMATCH);
			assertEquals(c0.param, -1);

			QualityScoreTreatment treatment = policy1.treatment;
			assertNotNull(treatment);

			assertThat(treatment.type, is(QualityScoreTreatmentType.PRESERVE));
			assertThat(treatment.param, is(40));
		}
	}

}

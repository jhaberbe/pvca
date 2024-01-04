import pymer4
import pathlib 
import pandas as pd
import scipy.linalg
from tqdm import tqdm

class PCVA:
    """Principal Variance Component Analysis
    """

    def __init__(self, expression_data: pd.DataFrame, experimental_data: pd.DataFrame):
        """Constructor

        Args:
            expression_data (pd.DataFrame): Expression data in (samples, genes) shape
            experimental_data (pd.DataFrame): Experimental data in (samples, features) shape
        """
        self.expression_data = expression_data
        self.experimental_data = experimental_data

        self._eigenvalues = None
        self._eigenvectors = None

        self.random_effects_matrix = None
        self.normalized_variance = None

        self._centered_expression_data = False

    def run(self, formula: str, sample_column: str = "sample", variance_explained_threshold: float = 0.90) -> pd.DataFrame:
        if not self._centered_expression_data:
            self.center_expression_data()
            self._centered_expression_data = True

        self.pca_expression()

        return self.run_linear_mixed_effects_model(formula, sample_column = sample_column, variance_explained_threshold = variance_explained_threshold)
    
    def center_expression_data(self) -> None:
        self.expression_data = (self.expression_data.T - self.expression_data.T.mean(axis=0)).T
    
    def pca_expression(self) -> None:
        self._eigenvalues, self._eigenvectors = scipy.linalg.eig(self.expression_data.corr())
    
    def run_linear_mixed_effects_model(self, formula: str, sample_column: str = "sample", variance_explained_threshold: float = 0.90) -> pd.DataFrame:
        assert formula.startswith("component"), ValueError("The formula must specify 'PC' as the dependent variable (ie. 'PC ~ 1 + feature1 + feature2 ...'")

        # We first select the number of components that explain at least our threshold of variance (minimum 90%)
        # scale the eigenvalues to sum to 1
        eigen_values_scaled = self._eigenvalues / self._eigenvalues.sum()
        # select the number of PCs to meet the threshold
        print(eigen_values_scaled.astype(float))
        print((eigen_values_scaled.cumsum().astype(float) < variance_explained_threshold))
        n_principal_components_included = (eigen_values_scaled.cumsum().astype(float) < variance_explained_threshold).sum() + 1


        # Construct the table we will be using to run the linear mixed effects model
        pc_data_matrix = pd.Series(self._eigenvectors[:, :n_principal_components_included].T.flatten()).astype(float)
        sample_components = pd.DataFrame(dict(
            sample = list(range(1, self.experimental_data.shape[0] + 1)) * n_principal_components_included,
            component = pc_data_matrix,
        ))
        sample_components = pd.merge(sample_components, self.experimental_data, how="left", on=sample_column)
        sample_components["PC"] = [i // self.experimental_data.shape[0] for i in range(sample_components.shape[0])]


        # Run the linear mixed effects model
        self.random_effects_matrix = pd.DataFrame()

        # for each principal component
        for pc_index in tqdm(range(n_principal_components_included)):
            # we model the experimental data to the value of the principal component
            pc_data_matrix = sample_components.query(f"PC == {pc_index}").copy()

            # we run the linear mixed effects model on the data to our PC value
            model = pymer4.models.Lmer(formula, data=pc_data_matrix)
            model.fit()

            # and collect the variance explained into this random_effects_matrix
            self.random_effects_matrix[pc_index] = model.ranef_var["Var"]

        # clean up the random effects matrix so that the variance explained is normalized
        self.normalized_variance = self.random_effects_matrix.T.apply(lambda x: x / x.sum(), axis=1)

        for pc_index in self.normalized_variance.index:
            self.normalized_variance.loc[pc_index] = eigen_values_scaled[pc_index] * normalized_variance.loc[pc_index]

        return self.normalized_variance.astype(float).sum()
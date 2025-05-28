# [Annotating cell clusters in single cell RNA-seq datasets](https://pluto.bio/resources/Learning%20Series/annotating-clusters-in-scrnaseq)

Single cell RNA sequencing (scRNA-seq) has revolutionized the way we study gene expression at the individual cell level. However, once you’ve performed clustering to group similar cells together, you’re faced with one of the most challenging tasks in scRNA-seq analysis: annotating your clusters. Cluster annotation is the process of assigning biological meaning to these groups, essentially identifying the cell types or states that each cluster represents.

For biologists who may not have much experience in bioinformatics, this task can feel daunting. But don’t worry! In this blog post, we’ll walk you through various strategies and tools that can help you confidently annotate your scRNA-seq clusters. From leveraging existing gene sets and reference databases to using advanced AI-powered tools, we’ve got you covered.

While there are many advanced algorithms and AI-powered tools that can aid in annotation, the most accurate and reliable approach often begins with a biology-first mindset. By combining your deep understanding of biology with available computational resources, you can enhance the reliability of your annotations and uncover novel insights that might otherwise be missed. Let’s dive into the key strategies and tools available for cluster annotation.

## [List of annotation tools and approaches](https://airtable.com/appMd0h4vP7gzQaeK/shrgmvY3ZvswENjkJ/tblgv3JRYlbD34DYD)
This list is just a starting place, there are many more tools available! Feel free to submit an approach to have it added to our growing list.

## Approaches for cell type annotation
1. Use existing gene sets or reference databases
One of the most straightforward ways to annotate your scRNA-seq clusters is to use existing gene sets or reference expression profiles. These are predefined collections of genes that are associated with specific cell types or biological processes.
  - Gene Ontology (GO) terms and Hallmark gene sets: These are well-established gene sets that describe common biological functions or pathways. By comparing the expression of genes in your clusters to these sets, you can gain insights into what each cluster might represent, such as identifying whether a cluster is related to immune response or cell cycle activity.
  - Reference datasets and cell type databases: For more precise cell type annotation, you can leverage curated reference datasets and cell type databases that contain gene expression profiles for specific cell types, tissues, and conditions. These datasets are invaluable for annotating scRNA-seq clusters by comparing their gene expression signatures to known profiles. Some of the most widely used resources include:
    - Human Cell Atlas (HCA): The Human Cell Atlas is a comprehensive reference dataset that provides gene expression data from a broad range of human cell types across multiple tissues. It's one of the most widely used reference resources for mapping and annotating human cell types.
    - PanglaoDB: PanglaoDB is a freely accessible database that compiles a list of marker genes for different human and mouse cell types, offering a useful reference for identifying clusters based on known cell markers.
    - Allen Brain Atlas: The Allen Brain Atlas provides detailed gene expression data specifically for neural and brain-related studies. It’s particularly valuable when working with clusters derived from brain tissues or studying neuronal subtypes.
    - MapMyCells: MapMyCells is a tool that allows you to map your scRNA-seq data to reference cell type atlases, such as the Human Cell Atlas. This enables you to match your clusters against well-characterized cell types from various tissues or conditions, providing fast and accurate annotations for your data.
    > Tip: These reference datasets are powerful, but they may not capture all the subtle biological differences or novel cell types present in your dataset. Combining these references with other methods, like AI-powered tools or manual annotation, can lead to more accurate results.

2. Leverage AI-powered annotation tools
As AI continues to make its mark in bioinformatics, several AI-based tools are now available to assist with cluster annotation. These tools are designed to predict cell types by analyzing gene expression patterns and comparing them to large-scale, annotated datasets
  - scGPT: One such tool is scGPT (Single-cell GPT), an AI-powered tool that uses generative pre-trained transformers (GPT) to predict cell types based on RNA-seq data. This AI model has been trained on vast datasets, enabling it to recognize patterns across a wide variety of cell types and tissues. With just a few clicks, scGPT can provide you with annotation suggestions for your clusters.
  - Pluto AI: Pluto offers an integrated AI-powered annotation tool that functions as a "copilot" to your biological domain expertise and offers initial ideas for cell types. A large language model (LLM) reads the marker genes for each cluster and suggests possible cell types (including a confidence of high/med/low) and provides rationale for why it is suggesting those cell types. For example, there may be some characteristic immune-related genes or other cell type-defining genes that are expressed in a cluster.

3. Apply machine learning-based approaches for cell type annotation
There are several machine learning-based tools that perform classification of cell types by comparing gene expression data to reference datasets. These tools typically use supervised learning models that have been trained on specific datasets and can offer more specific and targeted predictions.
  - CellTypist: This tool uses a supervised machine learning approach to classify cells into types based on expression profiles. It is trained on large-scale single-cell RNA-seq datasets and can predict cell types for a variety of tissues and conditions. You can also train CellTypist on your own data or use its existing models for annotations.
  - scMatch: Similar to CellTypist, scMatch uses a machine learning classifier to annotate single-cell clusters based on predefined cell type labels. It compares your dataset to a reference set of known cell types and makes predictions based on the closest matches. It's particularly useful for studies where you have a smaller set of reference cell types or marker genes.
  > Tip: Both AI and machine learning tools can provide rapid, data-driven annotations. However, it's important to cross-check the predictions, especially in cases where your data may contain novel or underrepresented cell types that might not be present in the reference models.

4. Use other algorithmic methods or specialized software packages
Another approach to cluster annotation involves using algorithmic methods and specialized bioinformatics software packages. These methods typically rely on statistical algorithms to identify the most likely cell type based on gene expression data.
  - scCatch: This tool uses a computational approach to predict cell types based on their gene expression profiles. scCatch is particularly useful if you want to annotate your clusters using a small set of known marker genes or if you want to integrate multiple datasets to get more robust annotations.
  - Seurat and SingleR: Popular bioinformatics packages like Seurat (for clustering analysis) and SingleR (for reference-based annotation) can also be extremely helpful. Seurat offers various ways to visualize and cluster your data, while SingleR performs annotation by comparing your clusters’ gene expression to reference cell type libraries.
  - Marker gene-based approaches: Many tools allow you to annotate clusters by looking at known marker genes. For example, SCINA uses a set of marker genes to classify cell types. You can also manually inspect the expression of these markers in your clusters to assign annotations.
  > Tip: For a more flexible approach, consider combining multiple algorithms to cross-check the results and ensure the accuracy of your annotations.

5. Apply domain expertise and manual annotation
While computational tools are invaluable for streamlining the annotation process, domain expertise and manual annotation should remain at the core of your analysis. Your deep understanding of the biological context of your experiment—whether it’s tissue-specific cell types, cell differentiation states, or pathology-related markers—cannot be replaced by automated tools.

## Why biology-first is the best annotation approach
- Contextual knowledge: Biological expertise allows you to consider the experimental context and biological relevance of your clusters. For example, if you're studying a disease state or developmental process, you’ll have a sense of what cell types or cell states should be present, and this knowledge can guide the annotation process.
- Identifying novel cell types: Automated methods are great for classifying known cell types, but manual annotation allows you to make informed decisions when faced with novel or uncharacterized cell types that might not be well-represented in reference databases or AI models.
- Identifying specific cell types: Automated methods are great for classifying broad cell types, but manual annotation allows for increased granularity to identify specific cell types. For example, automated analysis of T cells may be easily achieved but classification of a subset of CD4+ memory T cells may require manual annotation.
- Biological validation: Using domain knowledge to validate the predictions made by AI tools or algorithms helps ensure that the annotations are biologically meaningful. After all, computational tools are often trained on reference datasets, which might not capture the full complexity or variability of your own dataset.

## How to approach manual annotation
- Examine marker genes: Use known cell type markers to check the expression of specific genes in your clusters. This allows you to validate the identity of the cells and make manual adjustments if necessary.
- Work with colleagues or collaborators: If you're unsure about certain clusters, consult with other experts in the field. Domain experts can provide insights into the biological significance of your clusters and guide you in making informed annotations.

## Other important considerations for manual annotation
- Take detailed notes about your rationale for selected annotation: Unlike applying a single tool, manual annotation is not end-to-end reproducible. Rather, you should record your rationale for annotating a particular cluster, including any tools or reference data sources that were used to arrive at that annotation.
- While manual annotation may take more time, it ensures that you capture the full biological complexity of your dataset.

> Tip: Manual annotation can be especially important in cases where you're studying novel or poorly characterized cell types that aren't well represented in existing gene sets or reference databases.

Combining strategies for best results
The best approach to annotating your scRNA-seq clusters often involves a combination of methods. Start with biological expertise and manual review to provide a solid foundation, then supplement with computational tools (like MapMyCells, AI-powered tools, or machine learning classifiers) for efficiency and to confirm your hypotheses. This layered approach allows you to combine biological validation with computational efficiency.

Conclusion: Making cell type annotation less daunting
At the end of the day, the most reliable way to annotate your scRNA-seq clusters is by combining biological expertise with computational tools. While AI and machine learning methods can help speed up the process and provide initial annotations, manual annotation driven by your understanding of biology is essential for ensuring that your results are meaningful and accurate.

By starting with a biology-first approach and using computational tools to complement your expertise, you’ll be able to make more confident, biologically accurate annotations that reveal true insights about your data. Remember, the goal of annotation is to uncover the biology within your dataset—so trust your knowledge and leverage the best tools to help you achieve that!

At Pluto, we strive to make the bioinformatics process as user-friendly as possible. Our platform offers easy-to-use tools for cluster annotation, leveraging the power of AI and established reference databases to support your research. We hope this guide helps you feel more confident as you annotate your clusters and move forward with your scRNA-seq analysis!

Next steps
Explore our annotation tool: Check out our built-in annotation tool for a simple, automated approach to cluster annotation.

Learn more: Dive deeper into scRNA-seq analysis with our tutorials and resources.

Have any questions or need help with annotating your scRNA-seq clusters? Reach out to our support team or check out our knowledge base for more resources. We’re here to help you every step of the way!

List of scRNA-seq annotation tools and approaches
We'd like to thank Matthew Bernstein (@Matthew_N_B) for sharing a fantastic list of single-cell RNA-seq annotation tools, which we've incorporated into the following table.
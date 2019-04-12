package edu.cg;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.util.function.Supplier;

public class SeamsCarver extends ImageProcessor {

	// MARK: An inner interface for functional programming.
	@FunctionalInterface
	interface ResizeOperation {
		BufferedImage resize();
	}

	// MARK: Fields
	private int numOfSeams;
	private ResizeOperation resizeOp;
	boolean[][] imageMask;
	int[][] indexMatrix;
	long[][] energyMatrix;
	int[][] grayScaleMatrix;
	int currWidth;
	boolean[][] outImageMask;

	public SeamsCarver(Logger logger, BufferedImage workingImage, int outWidth, RGBWeights rgbWeights,
			boolean[][] imageMask) {
		super((s) -> logger.log("Seam carving: " + s), workingImage, rgbWeights, outWidth, workingImage.getHeight());

		numOfSeams = Math.abs(outWidth - inWidth);
		this.imageMask = imageMask;
		if (inWidth < 2 | inHeight < 2)
			throw new RuntimeException("Can not apply seam carving: workingImage is too small");

		if (numOfSeams > inWidth / 2)
			throw new RuntimeException("Can not apply seam carving: too many seams...");

		// Setting resizeOp by with the appropriate method reference
		if (outWidth > inWidth)
			resizeOp = this::increaseImageWidth;
		else if (outWidth < inWidth)
			resizeOp = this::reduceImageWidth;
		else
			resizeOp = this::duplicateWorkingImage;

		indexMatrix = createIndexMatrix();
		grayScaleMatrix = createGrayScaleValuesMatrix();
		outImageMask = new boolean[inHeight][outWidth];
		createEnergyMatrix();

		this.logger.log("preliminary calculations were ended.");
	}

	public BufferedImage resize() {
		return resizeOp.resize();
	}

	private BufferedImage reduceImageWidth() {
		return reduceImage(inWidth - outWidth, this::createNewImageFromIndexMatrix);
	}

	private BufferedImage increaseImageWidth() {
		return reduceImage(outWidth - inWidth, this::buildImage);
	}

	private BufferedImage reduceImage(int numberOfSeamsToRemove, Supplier<BufferedImage> image) {
		for (int i = 0; i < numberOfSeamsToRemove; i++) {
			int[] seam = findSeamToRemove();
			removeSeamFromIndexMatrix(seam);
			removeSeamFromEnergyMatrix(seam);
			currWidth--;
		}

		return image.get();
	}

	public BufferedImage showSeams(int seamColorRGB) {
		int newWidthSize = Math.abs(inWidth - outWidth);
		BufferedImage output = duplicateWorkingImage();

		for (int i = 0; i < newWidthSize; i++) {
			int[] seam = findSeamToRemove();
			drawSeamOnImage(seamColorRGB, seam, output);
			removeSeamFromIndexMatrix(seam);
			removeSeamFromEnergyMatrix(seam);
			currWidth--;
		}

		return output;
	}

	public boolean[][] getMaskAfterSeamCarving() {
		if (outWidth == inWidth) {
			return imageMask;
		}

		if (outWidth < inWidth) {
			setForEachParameters(outWidth, inHeight);

			forEach((y, x) -> {
				outImageMask[y][x] = imageMask[y][x];
			});
		}

		// inWidth > outWidth was handled in "buildImage".

		return outImageMask;
	}

	private int[][] createIndexMatrix() {
		int[][] matrix = new int[inHeight][inWidth];
		currWidth = inWidth;
		setForEachInputParameters();

		forEach((y, x) -> {
			matrix[y][x] = x;
		});

		return matrix;
	}

	private int[][] createGrayScaleValuesMatrix() {
		BufferedImage grayScaledImage = greyscale();
		int[][] grayScaleMatrix = new int[inHeight][inWidth];
		setForEachInputParameters();

		forEach((y, x) -> {
			grayScaleMatrix[y][x] = new Color(grayScaledImage.getRGB(x, y)).getBlue();
		});

		return grayScaleMatrix;
	}

	private void createEnergyMatrix() {
		energyMatrix = new long[inHeight][inWidth];
		setForEachInputParameters();

		forEach(this::calculatePixelsEnergy);
	}

	private int[] findSeamToRemove() {
		long[][] costMatrix = createCostMatrix();
		int[] seam = new int[inHeight];

		int minIndex = 0;
		long minValue = Long.MAX_VALUE;

		// Finding the minimum value in the last row.
		for (int j = 0; j < currWidth; j++) {
			if (costMatrix[inHeight - 1][indexMatrix[inHeight - 1][j]] < minValue) {
				minIndex = j;
				minValue = costMatrix[inHeight - 1][indexMatrix[inHeight - 1][j]];
			}
		}

		seam[inHeight - 1] = minIndex;
		int j = minIndex;

		// Backtracking to find the seam.
		for (int i = inHeight - 1; i > 0; i--) {
			// j
			if (costMatrix[i][indexMatrix[i][j]] == energyMatrix[i][indexMatrix[i][j]] + calcCU(i, j, costMatrix)) {
				seam[i - 1] = j;

				continue;
			}

			// j-1
			if (j != 0 && costMatrix[i][indexMatrix[i][j]] == energyMatrix[i][indexMatrix[i][j]]
					+ calcCL(i, j, costMatrix)) {
				seam[i - 1] = --j;

				continue;
			}

			// j+1
			seam[i - 1] = ++j;
		}

		return seam;
	}

	private long[][] createCostMatrix() {
		long[][] costMatrix = new long[inHeight][inWidth];
		setForEachParameters(currWidth, inHeight);

		forEach((y, x) -> {
			long CL, CU, CR;

			// If on the first row
			if (y == 0) {
				costMatrix[y][indexMatrix[y][x]] = energyMatrix[y][indexMatrix[y][x]];
				return;
			}

			CU = calcCU(y, x, costMatrix);

			// If on the first column
			if (x == 0) {
				CR = calcCR(y, x, costMatrix);

				costMatrix[y][indexMatrix[y][x]] = energyMatrix[y][indexMatrix[y][x]] + Math.min(CU, CR);
				return;
			}

			// If on the last column
			if (x == currWidth - 1) {
				CL = calcCL(y, x, costMatrix);

				costMatrix[y][indexMatrix[y][x]] = energyMatrix[y][indexMatrix[y][x]] + Math.min(CU, CL);
				return;
			}

			CL = calcCL(y, x, costMatrix);
			CR = calcCR(y, x, costMatrix);

			costMatrix[y][indexMatrix[y][x]] = energyMatrix[y][indexMatrix[y][x]] + Math.min(CL, Math.min(CU, CR));
		});

		return costMatrix;
	}

	private long calcCU(int y, int x, long[][] costMatrix) {
		return costMatrix[y - 1][indexMatrix[y - 1][x]] + getMiddleEdgeEnergy(y, x);
	}

	private long calcCR(int y, int x, long[][] costMatrix) {
		return costMatrix[y - 1][indexMatrix[y - 1][x + 1]]
				+ Math.abs(energyMatrix[y - 1][indexMatrix[y - 1][x]] - energyMatrix[y][indexMatrix[y][x + 1]])
				+ getMiddleEdgeEnergy(y, x);
	}

	private long calcCL(int y, int x, long[][] costMatrix) {
		return costMatrix[y - 1][indexMatrix[y - 1][x - 1]]
				+ Math.abs(energyMatrix[y - 1][indexMatrix[y - 1][x]] - energyMatrix[y][indexMatrix[y][x - 1]])
				+ getMiddleEdgeEnergy(y, x);
	}

	private long getMiddleEdgeEnergy(int y, int x) {
		if (x == 0 || x == currWidth - 1) {
			return 0;
		}

		return Math.abs(energyMatrix[y][indexMatrix[y][x - 1]] - energyMatrix[y][indexMatrix[y][x + 1]]);
	}

	private void removeSeamFromEnergyMatrix(int[] seam) {
		setForEachParameters(currWidth, inHeight);
		int j;

		for (int i = 0; i < inHeight; i++) {
			j = seam[i];

			if (i > 0) {
				calculatePixelsEnergy(i - 1, j);
			}

			if (j > 0 && j < currWidth - 1) {
				calculatePixelsEnergy(i, j - 1);
			}
		}
	}

	private void calculatePixelsEnergy(int y, int x) {
		int currPixelValue = grayScaleMatrix[y][indexMatrix[y][x]];

		int e1 = Math.abs(currPixelValue - (y == (inHeight - 1) ? grayScaleMatrix[y - 1][indexMatrix[y - 1][x]]
				: grayScaleMatrix[y + 1][indexMatrix[y + 1][x]]));
		int e2 = Math.abs(currPixelValue - (x == (currWidth - 1) ? grayScaleMatrix[y][indexMatrix[y][x - 1]]
				: grayScaleMatrix[y][indexMatrix[y][x + 1]]));
		long e3 = imageMask[y][indexMatrix[y][x]] ? Integer.MAX_VALUE : 0;

		energyMatrix[y][indexMatrix[y][x]] = e1 + e2 + e3;
	}

	private BufferedImage createNewImageFromIndexMatrix() {
		BufferedImage output = newEmptyOutputSizedImage();
		setForEachOutputParameters();

		forEach((y, x) -> {
			output.setRGB(x, y, workingImage.getRGB(indexMatrix[y][x], y));
		});

		return output;
	}

	private void removeSeamFromIndexMatrix(int[] seam) {
		for (int i = 0; i < inHeight; i++) {
			if(seam[i] + 1 == currWidth) {
				indexMatrix[i][seam[i]] = 0;
			}
			
			for (int j = seam[i] + 1; j < currWidth; j++) {
				indexMatrix[i][j - 1] = indexMatrix[i][j];
				imageMask[i][j - 1] = imageMask[i][j];
			}
		}
	}

	private BufferedImage buildImage() {
		BufferedImage ans = newEmptyOutputSizedImage();
		int indexForMapping;

		for (int i = 0; i < inHeight; i++) {
			indexForMapping = 0;

			for (int j = 0; j < inWidth; j++) {
				int pixelColor = workingImage.getRGB(j, i);
				boolean isMasked = imageMask[i][j];

				ans.setRGB(j + indexForMapping, i, pixelColor);
				outImageMask[i][j + indexForMapping] = isMasked;

				// If the pixel is not in the index matrix
				// the pixel was remove with a seam => need to duplicate.
				if (j != indexMatrix[i][j - indexForMapping]) {
					ans.setRGB(j + (++indexForMapping), i, pixelColor);
					outImageMask[i][j + indexForMapping] = isMasked;
				}
			}
		}

		return ans;
	}

	private void drawSeamOnImage(int seamColorRGB, int[] seam, BufferedImage output) {
		for (int i = 0; i < inHeight; i++) {
			output.setRGB(indexMatrix[i][seam[i]], i, seamColorRGB);
		}
	}
}

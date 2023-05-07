using System;
using UnityEditor;
using UnityEngine;
using UnityEngine.XR;
using Vectors;

[ExecuteAlways]
[RequireComponent(typeof(VectorRenderer))]
public class ModelScript : MonoBehaviour {
    
    [NonSerialized] 
    private VectorRenderer vectors;

    [SerializeField, Range(0f, 1f)] 
    private float Time;
    
    [SerializeField,HideInInspector]internal Matrix4x4 matrixA = Matrix4x4.identity;
    [SerializeField,HideInInspector]internal Matrix4x4 matrixB = Matrix4x4.identity;
    [SerializeField,HideInInspector]internal Matrix4x4 matrixI = Matrix4x4.identity;

    [SerializeField] 
    private bool interpolateTranslation = true;
    [SerializeField] 
    private bool interpolateRotation = true;
    [SerializeField] 
    private bool interpolateScale = true;

    void OnEnable() {
        vectors = GetComponent<VectorRenderer>();
    }

    void Update()
    {

        if (interpolateTranslation) // Check if user wants to interpolate translation
        {
            var interpolatedTranslation = (1.0f - Time) * GetTranslationVector(matrixA) + Time * GetTranslationVector(matrixB); // Interpolate translation
            SetTranslationVector(ref matrixI, interpolatedTranslation); // Set the translation of interpolated matrix to interpolated translation
        }
        else
        {
            SetTranslationVector(ref matrixI, GetTranslationVector(matrixA)); // If translation is ignored,
                                                                              // set translation of interpolated matrix to starting matrix
        }

        if (interpolateRotation) // Check if user wants to interpolate rotation
        {
            var rotationA = GetQuaternionFromMatrix(matrixA); // Extract quaternion from matrix A
            var rotationB = GetQuaternionFromMatrix(matrixB); // Extract quaternion from matrix B
            
            var rotationI = GetInterpolatedQuaternion(rotationA, rotationB); // Use quaternion A and quaternion B as
                                                                                      // parameters for the interpolate quaternion function
                                                                                      
            var rotationMatrixI = GetMatrixFromQuaternion(rotationI); // Convert the interpolated quaternion into a rotation matrix
            matrixI = SetRotation( matrixI, rotationMatrixI); // Use the rotation matrix as the parameter for the set rotation function
        }
        else
        {
            matrixI = SetRotation(matrixI, Matrix4x4.identity); // If rotation is ignored, use the identity matrix
                                                                                // as parameter for the set rotation function
        }

        if (interpolateScale) // Check if user wants to interpolate scale
        {
            SetScale(ref matrixI, GetInterpolatedScale(matrixA, matrixB)); // Use the return of the interpolated scale magnitudes
                                                                                     // as parameter for the set scale function
        }
        else
        {
            SetScale(ref matrixI, new Vector3(1,1,1)); // If scale is ignored, use a vector with all components
                                                            // as ones as parameter for the set scale function
        }

        // Variables for vector addition used for the visualization of the base vectors per matrix
        
        Vector3 AX = matrixA.GetColumn(0); // GetColumn returns a vector representing the column with index, example index:0 returns (m00, m10, m20)
        Vector3 AY = matrixA.GetColumn(1); 
        Vector3 AZ = matrixA.GetColumn(2);
        
        Vector3 BX = matrixB.GetColumn(0);
        Vector3 BY = matrixB.GetColumn(1);
        Vector3 BZ = matrixB.GetColumn(2);
        
        Vector3 IX = matrixI.GetColumn(0);
        Vector3 IY = matrixI.GetColumn(1);
        Vector3 IZ = matrixI.GetColumn(2);
        
        // The matrices translation vectors
        Vector3 TA = GetTranslationVector(matrixA);
        Vector3 TB = GetTranslationVector(matrixB);
        Vector3 TI = GetTranslationVector(matrixI);
        
        using (vectors.Begin()) {
            
            //Visualize matrix A
            vectors.Draw(TA, TA+AX, Color.red);
            vectors.Draw(TA, TA+AY, Color.green);
            vectors.Draw(TA, TA+AZ, Color.blue);
            
            //Visualize matrix B
            vectors.Draw(TB, TB+BX, Color.red);
            vectors.Draw(TB, TB+BY, Color.green);
            vectors.Draw(TB, TB+BZ, Color.blue);
            
            //Visualize box from interpolated matrix
            vectors.Draw(TI, TI+IX, Color.red);
            vectors.Draw(TI+IY, TI+IX+IY, Color.red);
            vectors.Draw(TI+IZ, TI+IX+IZ, Color.red);
            vectors.Draw(TI+IY+IZ, TI+IX+IY+IZ, Color.red);
            
            vectors.Draw(TI, TI+IY, Color.green);
            vectors.Draw(TI+IX, TI+IY+IX, Color.green);
            vectors.Draw(TI+IZ, TI+IY+IZ, Color.green);
            vectors.Draw(TI+IX+IZ, TI+IY+IX+IZ, Color.green);
            
            vectors.Draw(TI, TI+IZ, Color.blue);
            vectors.Draw(TI+IX, TI+IZ+IX, Color.blue);
            vectors.Draw(TI+IY, TI+IZ+IY, Color.blue);
            vectors.Draw(TI+IX+IY, TI+IZ+IX+IY, Color.blue);
        }
    }

    public float GetDeterminantFromMatrix(Matrix4x4 matrix)
    {
        // Determinant by expanding the cofactors from the recursive definition of a 4x4 matrix determinant
        var determinant = matrix.m00 * ((matrix.m11 * (matrix.m22 * matrix.m33 - matrix.m23 * matrix.m32)) + 
                                    (matrix.m12 * (matrix.m23 * matrix.m31 - matrix.m21 * matrix.m33)) +
                                    (matrix.m13 * (matrix.m21 * matrix.m32 - matrix.m22 * matrix.m31))) -
                      matrix.m01 * ((matrix.m10 * (matrix.m22 * matrix.m33 - matrix.m23 * matrix.m32)) +
                                    (matrix.m12 * (matrix.m23 * matrix.m30 - matrix.m20 * matrix.m33)) +
                                    (matrix.m13 * (matrix.m20 * matrix.m32 - matrix.m22 * matrix.m30))) +
                      matrix.m02 * ((matrix.m10 * (matrix.m21 * matrix.m33 - matrix.m23 * matrix.m31)) +
                                    (matrix.m11 * (matrix.m23 * matrix.m30 - matrix.m20 * matrix.m33)) +
                                    (matrix.m13 * (matrix.m20 * matrix.m31 - matrix.m21 * matrix.m30))) -
                      matrix.m03 * ((matrix.m10 * (matrix.m21 * matrix.m32 - matrix.m22 * matrix.m31)) +
                                    (matrix.m11 * (matrix.m22 * matrix.m30 - matrix.m20 * matrix.m32)) +
                                    (matrix.m12 * (matrix.m20 * matrix.m31 - matrix.m21 * matrix.m30)));
        return determinant;
    }

    public Vector3 GetTranslationVector(Matrix4x4 matrix)
    {
        return new Vector3(matrix.m03, matrix.m13, matrix.m23); // Same as matrix.GetColumn(3)
    }
    
    public void SetTranslationVector(ref Matrix4x4 matrix, Vector3 translationVector)
    {
        matrix.m03 = translationVector.x; // Same as matrix.SetColumn(3, translationVector)
        matrix.m13 = translationVector.y;
        matrix.m23 = translationVector.z;
    }

    public Vector3 GetInterpolatedScale(Matrix4x4 matrixA, Matrix4x4 matrixB)
    {
        // Get magnitude of column vectors for matrix A
        var scaleAX = MathF.Sqrt(matrixA.m00 * matrixA.m00 + matrixA.m10 * matrixA.m10 + matrixA.m20 * matrixA.m20);
        var scaleAY = MathF.Sqrt(matrixA.m01 * matrixA.m01 + matrixA.m11 * matrixA.m11 + matrixA.m21 * matrixA.m21);
        var scaleAZ = MathF.Sqrt(matrixA.m02 * matrixA.m02 + matrixA.m12 * matrixA.m12 + matrixA.m22 * matrixA.m22);
        
        // Get magnitude of column vectors for matrix B
        var scaleBX = MathF.Sqrt(matrixB.m00 * matrixB.m00 + matrixB.m10 * matrixB.m10 + matrixB.m20 * matrixB.m20);
        var scaleBY = MathF.Sqrt(matrixB.m01 * matrixB.m01 + matrixB.m11 * matrixB.m11 + matrixB.m21 * matrixB.m21);
        var scaleBZ = MathF.Sqrt(matrixB.m02 * matrixB.m02 + matrixB.m12 * matrixB.m12 + matrixB.m22 * matrixB.m22);
        
        // Interpolate scale magnitudes
        var scaleIX = (1.0f - Time) * scaleAX + Time * scaleBX;
        var scaleIY = (1.0f - Time) * scaleAY + Time * scaleBY;
        var scaleIZ = (1.0f - Time) * scaleAZ + Time * scaleBZ;
        
        return new Vector3(scaleIX, scaleIY, scaleIZ); //Scaling vector
    }

    public Vector3 GetScaleFromMatrix(Matrix4x4 matrix)
    {
        // Get magnitude of column vectors for matrix
        var scaleX = MathF.Sqrt(matrix.m00 * matrix.m00 + matrix.m10 * matrix.m10 + matrix.m20 * matrix.m20);
        var scaleY = MathF.Sqrt(matrix.m01 * matrix.m01 + matrix.m11 * matrix.m11 + matrix.m21 * matrix.m21);
        var scaleZ = MathF.Sqrt(matrix.m02 * matrix.m02 + matrix.m12 * matrix.m12 + matrix.m22 * matrix.m22);
        
        return new Vector3(scaleX, scaleY, scaleZ); //Scaling vector
    }

    public void SetScale(ref Matrix4x4 matrix, Vector3 scaleVector)
    {
        matrix.SetColumn(0,matrix.GetColumn(0)*scaleVector.x); // Multiply matrix's X column with scale-vector's X component
        matrix.SetColumn(1,matrix.GetColumn(1)*scaleVector.y); // Multiply matrix's Y column with scale-vector's Y component
        matrix.SetColumn(2,matrix.GetColumn(2)*scaleVector.z); // Multiply matrix's Z column with scale-vector's Z component
    }

    public Quaternion GetQuaternionFromMatrix(Matrix4x4 matrix)
    {
        var matrixColumnX = matrix.GetColumn(0);
        var matrixColumnY = matrix.GetColumn(1);
        var matrixColumnZ = matrix.GetColumn(2);
        
        // Normalize column vectors of matrix
        matrix.SetColumn(0, GetNormalizedVector(matrixColumnX));
        matrix.SetColumn(1, GetNormalizedVector(matrixColumnY));
        matrix.SetColumn(2, GetNormalizedVector(matrixColumnZ));
        
        // Components of output quaternion
        float w = 0, x = 0, y = 0, z = 0;
        
        // Find which is larger
        float wSubOne = matrix.m00 + matrix.m11 + matrix.m22;
        float xSubOne = matrix.m00 - matrix.m11 - matrix.m22;
        float ySubOne = matrix.m11 - matrix.m00 - matrix.m22;
        float zSubOne = matrix.m22 - matrix.m00 - matrix.m11;
        
        int biggestIndex = 0;
        float biggestSubOne = wSubOne;
        if (xSubOne > biggestSubOne)
        {
            biggestSubOne = xSubOne;
            biggestIndex = 1;
        }

        if (ySubOne > biggestSubOne)
        {
            biggestSubOne = ySubOne;
            biggestIndex = 2;
        }

        if (zSubOne > biggestSubOne)
        {
            biggestSubOne = zSubOne;
            biggestIndex = 3;
        }
        // Calculate the square root and divide
        float biggestVal = MathF.Sqrt(biggestSubOne + 1.0f) * 0.5f;
        float mult = 0.25f / biggestVal;
        // Compute the components based on which is larger
        switch (biggestIndex)
        {
            case 0:
                w = biggestVal;
                x = (matrix.m21 - matrix.m12) * mult;
                y = (matrix.m02 - matrix.m20) * mult;
                z = (matrix.m10 - matrix.m01) * mult;
                break;

            case 1:
                x = biggestVal;
                w = (matrix.m21 - matrix.m12) * mult;
                y = (matrix.m10 + matrix.m01) * mult;
                z = (matrix.m02 + matrix.m20) * mult;
                break;

            case 2:
                y = biggestVal;
                w = (matrix.m02 - matrix.m20) * mult;
                x = (matrix.m10 + matrix.m01) * mult;
                z = (matrix.m21 + matrix.m12) * mult;
                break;

            case 3:
                z = biggestVal;
                w = (matrix.m10 - matrix.m01) * mult;
                x = (matrix.m02 + matrix.m20) * mult;
                y = (matrix.m21 + matrix.m12) * mult;
                break;

        }
        
        return new Quaternion(x, y, z, w);
    }
    
    public Quaternion GetInterpolatedQuaternion(Quaternion quatA, Quaternion quatB)
    {
        // Variables for components of quatA and quatB
        float ax = quatA.x;
        float ay = quatA.y;
        float az = quatA.z;
        float aw = quatA.w;
        
        float bx = quatB.x;
        float by = quatB.y;
        float bz = quatB.z;
        float bw = quatB.w;

        float t = Time; // Interpolation variable

        float cx,cy,cz,cw; // Output components of interpolated quaternion
        
        float cosOmega = aw*bw + ax*bx + ay*by + az*bz; // Get the cosine of the angle between the quaternions
        
        if (cosOmega < 0.0f) { // If the cosine is negative, negate components of one of the quaternions to get shortest angle
            bw = -bw;
            bx = -bx;
            by = -by;
            bz = -bz;
            cosOmega = -cosOmega;
        }
        float k0, k1;
        if (cosOmega > 0.9999f) { // Check if close to avoid division by zero, to then use linear interpolation

            k0 = 1.0f-t;
            k1 = t;

        } else {

            float sinOmega = Mathf.Sqrt(1.0f - cosOmega*cosOmega); // Get the sin of the angle using the trigonometric one
            
            float omega = Mathf.Atan2(sinOmega, cosOmega); // Get the angle using sine and cosine in the Atan2 function

            float oneOverSinOmega = 1.0f / sinOmega; // Get inverse so division is only required one time

            // Get interpolation parameters
            k0 = Mathf.Sin((1.0f - t) * omega) * oneOverSinOmega;
            k1 = Mathf.Sin(t * omega) * oneOverSinOmega;
        }
        
        // Interpolation of the components
        cx = ax*k0 + bx*k1;
        cy = ay*k0 + by*k1;
        cz = az*k0 + bz*k1;
        cw = aw*k0 + bw*k1;

        return new Quaternion(cx, cy, cz, cw);
    }

    public Matrix4x4 GetMatrixFromQuaternion(Quaternion quat)
    {
        var newMatrix = Matrix4x4.identity;
        
        // Formula for conversion

        // Column 1
        newMatrix.m00 = 1 - 2 * (quat.y * quat.y) - 2 * (quat.z * quat.z);
        newMatrix.m01 = 2 * (quat.x * quat.y) - 2 * (quat.w * quat.z);
        newMatrix.m02 = 2 * (quat.x * quat.z) + 2 * (quat.w * quat.y);
        // Column 2
        newMatrix.m10 = 2 * (quat.x * quat.y) + 2 * (quat.w * quat.z);
        newMatrix.m11 = 1 - 2 * (quat.x * quat.x) - 2 * (quat.z * quat.z);
        newMatrix.m12 = 2 * (quat.y * quat.z) - 2 * (quat.w * quat.x);
        // Column 3
        newMatrix.m20 = 2 * (quat.x * quat.z) - 2 * (quat.w * quat.y);
        newMatrix.m21 = 2 * (quat.y * quat.z) + 2 * (quat.w * quat.x);
        newMatrix.m22 = 1 - 2 * (quat.x * quat.x) - 2 * (quat.y * quat.y);
        
        return newMatrix;
    }

    public Matrix4x4 SetRotation(Matrix4x4 matrix, Matrix4x4 rotationMatrix)
    {
        // Sets the upper left 3x3 of the matrix parameter to the upper left 3x3 of the rotation matrix parameter
        matrix.m00 = rotationMatrix.m00;
        matrix.m01 = rotationMatrix.m01;
        matrix.m02 = rotationMatrix.m02;
        matrix.m10 = rotationMatrix.m10;
        matrix.m11 = rotationMatrix.m11;
        matrix.m12 = rotationMatrix.m12;
        matrix.m20 = rotationMatrix.m20;
        matrix.m21 = rotationMatrix.m21;
        matrix.m22 = rotationMatrix.m22;

        return matrix;
    }

    public Vector3 GetNormalizedVector(Vector3 vector)
    {
        var magnitude = MathF.Sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
        return vector / magnitude;
    }
}

[CustomEditor(typeof(ModelScript))]
public class ModelGUI : Editor {
    void OnSceneGUI() {
        var ex = target as ModelScript;
        if (ex == null) return;

        EditorGUI.BeginChangeCheck();
        var newTranslateA = (Vector4)Handles.PositionHandle(ex.GetTranslationVector(ex.matrixA), ex.GetQuaternionFromMatrix(ex.matrixA)); // Handle for position
        var newScaleA = Handles.ScaleHandle(ex.GetScaleFromMatrix(ex.matrixA), ex.GetTranslationVector(ex.matrixA), // Handle for scaling
            ex.GetQuaternionFromMatrix(ex.matrixA));
        var newRotationA =
            Handles.RotationHandle(ex.GetQuaternionFromMatrix(ex.matrixA), ex.GetTranslationVector(ex.matrixA)); // Handle for rotating
        
        var newTranslateB = (Vector4)Handles.PositionHandle(ex.GetTranslationVector(ex.matrixB), ex.GetQuaternionFromMatrix(ex.matrixB)); // Handle for position
        var newScaleB = Handles.ScaleHandle(ex.GetScaleFromMatrix(ex.matrixB), ex.GetTranslationVector(ex.matrixB), // Handle for scaling
            ex.GetQuaternionFromMatrix(ex.matrixB), 0.5f);
        var newRotationB =
            Handles.RotationHandle(ex.GetQuaternionFromMatrix(ex.matrixB), ex.GetTranslationVector(ex.matrixB)); // Handle for rotating
        
        if (EditorGUI.EndChangeCheck())
        {
            var a = Matrix4x4.identity;
            var b = Matrix4x4.identity;
            Undo.RecordObject(target, "Vector Positions");
            newTranslateA.w = 1; // Fix w component
            newTranslateB.w = 1; // Fix w component

            a = new Matrix4x4(
                newRotationA * new Vector3(1,0,0) * newScaleA.x,
                newRotationA * new Vector3(0,1,0) * newScaleA.y,
                newRotationA * new Vector3(0,0,1) * newScaleA.z, newTranslateA); // Construct new matrix A based on handles
            b = new Matrix4x4(
                newRotationB * new Vector3(1,0,0) * newScaleB.x,
                newRotationB * new Vector3(0,1,0) * newScaleB.y,
                newRotationB * new Vector3(0,0,1) * newScaleB.z, newTranslateB); // Construct new matrix B based on handles

            ex.matrixA = a;
            ex.matrixB = b;

            EditorUtility.SetDirty(target);
        }
    }
    public override void OnInspectorGUI()
    {
        base.OnInspectorGUI();
        
        var ex = target as ModelScript;
        if (ex == null) return;

        EditorGUI.BeginChangeCheck();

        EditorGUILayout.BeginHorizontal();
        EditorGUILayout.PrefixLabel("Matrix A");
        EditorGUILayout.BeginVertical();

        var resultA = Matrix4x4.identity;
        for (var i = 0; i < 4; i++)
        {
            EditorGUILayout.BeginHorizontal();
            for (var j = 0; j < 4; j++)
            {
                resultA[i, j] = EditorGUILayout.FloatField((ex.matrixA[i, j]));
            }
            EditorGUILayout.EndHorizontal();
        }
        
        EditorGUILayout.EndVertical();
        EditorGUILayout.EndHorizontal();
        
        EditorGUILayout.Space();
        
        EditorGUILayout.BeginHorizontal();
        EditorGUILayout.PrefixLabel("Matrix B");
        EditorGUILayout.BeginVertical();

        var resultB = Matrix4x4.identity;
        for (var k = 0; k < 4; k++)
        {
            EditorGUILayout.BeginHorizontal();
            for (var l = 0; l < 4; l++)
            {
                resultB[k, l] = EditorGUILayout.FloatField((ex.matrixB[k, l]));
            }
            EditorGUILayout.EndHorizontal();
        }
        
        EditorGUILayout.EndVertical();
        EditorGUILayout.EndHorizontal();
        
        EditorGUILayout.Space();
        
        EditorGUILayout.BeginHorizontal();
        EditorGUILayout.PrefixLabel("Matrix I");
        EditorGUILayout.BeginVertical();

        var resultI = Matrix4x4.identity;
        for (var o = 0; o < 4; o++)
        {
            EditorGUILayout.BeginHorizontal();
            for (var p = 0; p < 4; p++)
            {
                resultI[o, p] = EditorGUILayout.FloatField((ex.matrixI[o, p]));
            }
            EditorGUILayout.EndHorizontal();
        }
        
        EditorGUILayout.EndVertical();
        EditorGUILayout.EndHorizontal();
        
        EditorGUILayout.Space();

        EditorGUILayout.BeginHorizontal();
        EditorGUILayout.PrefixLabel("Determinant A");
        EditorGUILayout.FloatField((ex.GetDeterminantFromMatrix(ex.matrixA)));
        EditorGUILayout.EndHorizontal();
        
        EditorGUILayout.Space();

        EditorGUILayout.BeginHorizontal();
        EditorGUILayout.PrefixLabel("Determinant B");
        EditorGUILayout.FloatField((ex.GetDeterminantFromMatrix(ex.matrixB)));
        EditorGUILayout.EndHorizontal();
        
        EditorGUILayout.Space();

        EditorGUILayout.BeginHorizontal();
        EditorGUILayout.PrefixLabel("Determinant I");
        EditorGUILayout.FloatField((ex.GetDeterminantFromMatrix(ex.matrixI)));
        EditorGUILayout.EndHorizontal();
        

        if (EditorGUI.EndChangeCheck())
        {
            Undo.RecordObject(ex, "Changed matrix");
            ex.matrixA = resultA;
            ex.matrixB = resultB;
            EditorUtility.SetDirty(ex);
        }
    }
}
